program main

    use mod_var
    use mod_fem
    use mod_function
    implicit none

    ! Mesh and FEM variables
    real(PR), allocatable :: nodes(:,:)
    integer, allocatable :: quad_connectivity(:,:)
    integer :: num_nodes, num_quad_elements
    integer :: num_dofs, num_edges_x, num_edges_y

    ! Mesh parameters
    integer :: nx, ny, nodes_x, nodes_y
    real(PR) :: Lx, Ly, dx, dy

    ! Sparse matrix
    type(triplet_matrix) :: K_triplet, M_triplet, A_system_triplet
    type(sparse_matrix) :: K_sparse, M_sparse, A_system_sparse

    ! Solution vectors
    real(PR), allocatable :: U_solution(:), F_vector(:), U_solution_old(:)
    real(PR), allocatable :: U_analytical(:),U_boundary_perfect_conductor(:), F_rhs(:)

    ! Boundary conditions
    integer, allocatable :: boundary_dofs(:)
    integer :: num_boundary_dofs

    ! Element arrays
    real(PR) :: Ke(4,4),Me(4,4), fe(4)
    real(PR) :: elem_coords(4,2)
    integer :: elem_nodes(4), elem_edges(4)

    ! Time-stepping parameters
    real(PR) :: dt, t_final, current_time
    integer :: time_step, num_time_steps

    ! Misc
    real(PR) :: start_time, end_time,tolerance,residual
    integer :: i, j, elem, ierr,max_iter
    integer :: ix, iy, node_id, elem_id
    logical, allocatable :: is_bnd(:)
    integer, dimension(4) :: s

    ! =================================================================
    ! 1. SETUP
    ! =================================================================
    nx = 150
    ny = 150
    Lx = 1.0_PR
    Ly = 1.0_PR
    tolerance=0.0000000001_PR
    max_iter=2000

    ! --- Time parameters ---
    t_final = 2.0_PR
    dt = 0.01_PR
    num_time_steps = nint(t_final / dt)

    nodes_x = nx + 1
    nodes_y = ny + 1
    num_nodes = nodes_x * nodes_y
    num_quad_elements = nx * ny
    num_edges_x = nx * (ny + 1)
    num_edges_y = (nx + 1) * ny
    num_dofs = num_edges_x + num_edges_y
    dx = Lx / real(nx, PR)
    dy = Ly / real(ny, PR)

    write(*,*) 'Static Patch Test with Nedelec Elements'
    write(*,*) '========================================='
    write(*,*) 'Mesh:', nx, 'x', ny, 'elements'
    write(*,*) 'Total DoFs (edges):', num_dofs

    allocate(nodes(num_nodes, 2))
    allocate(quad_connectivity(num_quad_elements, 4))

    ! Generate nodes and connectivity (same as before)
    node_id = 0
    do iy = 1, nodes_y
        do ix = 1, nodes_x
            node_id = node_id + 1
            nodes(node_id, 1) = real(ix-1, PR) * dx
            nodes(node_id, 2) = real(iy-1, PR) * dy
        end do
    end do
    elem_id = 0
    do iy = 1, ny
        do ix = 1, nx
            elem_id = elem_id + 1
            quad_connectivity(elem_id, 1) = (iy-1) * nodes_x + ix
            quad_connectivity(elem_id, 2) = (iy-1) * nodes_x + ix + 1
            quad_connectivity(elem_id, 3) = iy * nodes_x + ix + 1
            quad_connectivity(elem_id, 4) = iy * nodes_x + ix
        end do
    end do

    ! =================================================================
    ! 2. ASSEMBLY
    ! =================================================================
    allocate(U_solution(num_dofs), F_vector(num_dofs), U_solution_old(num_dofs), F_rhs(num_dofs))
    F_vector = 0.0_PR
    U_solution = 0.0_PR ! Initial condition A(t=0) = 0
    U_solution_old = 0.0_PR


    call triplet_init(K_triplet, num_dofs, num_dofs * 9)
    call triplet_init(M_triplet, num_dofs, num_dofs * 9)
    write(*,*) 'Starting sparse assembly...'
    call cpu_time(start_time)

    do elem = 1, num_quad_elements
        elem_nodes = quad_connectivity(elem, :)
        do i = 1, 4
            elem_coords(i, :) = nodes(elem_nodes(i), :)
        end do

        ! Edge numbering (this logic is correct)
        ix = mod(elem-1, nx) + 1
        iy = (elem-1) / nx + 1
        elem_edges(1) = (iy-1)*nx + ix
        elem_edges(2) = num_edges_x + ix*ny + iy
        elem_edges(3) = iy*nx + ix
        elem_edges(4) = num_edges_x + (ix-1)*ny + iy

        ! Call the element routine
        call process_nedelec_element(elem_nodes, elem_edges, elem_coords,Me, Ke, fe, ierr)
        if (ierr /= 0) stop 'Error in element processing!'

        ! Assemble K and F with sign correction

        s = [1, 1, -1, -1] ! signs for [bottom, right, top, left]
        do i = 1, 4
            F_vector(elem_edges(i)) = F_vector(elem_edges(i)) + s(i) * fe(i)
            do j = 1, 4
                call triplet_add_entry(K_triplet, elem_edges(i), elem_edges(j), s(i) * s(j) * Ke(i,j))
                call triplet_add_entry(M_triplet, elem_edges(i), elem_edges(j), s(i) * s(j) * Me(i,j))
            end do
        end do
    end do

    call triplet_to_csr(K_triplet, K_sparse)
    call triplet_to_csr(M_triplet, M_sparse)
    call triplet_free(K_triplet)
    call triplet_free(M_triplet)

    ! --- Form the system matrix for Backward Euler: A = K + (1/dt) * M ---
    call triplet_init(A_system_triplet, num_dofs, K_sparse%nnz + M_sparse%nnz)
    call add_scaled_csr_to_triplet(A_system_triplet, K_sparse, 1.0_PR)
    call add_scaled_csr_to_triplet(A_system_triplet, M_sparse, 1.0_PR / dt)
    call triplet_to_csr(A_system_triplet, A_system_sparse)
    call triplet_free(A_system_triplet)

    call cpu_time(end_time)
    write(*,'(A,F8.3,A)') 'Assembly time: ', end_time - start_time, ' seconds'

    ! =================================================================
    ! 3. DEFINE BOUNDARY CONDITIONS
    ! =================================================================
    allocate(U_analytical(num_dofs)) ! Still needed for error calculation if desired
    allocate(U_boundary_perfect_conductor(num_dofs))
    U_boundary_perfect_conductor = 0.0_PR

    allocate(boundary_dofs(2*nx + 2*ny))
    num_boundary_dofs = 0
    do ix = 1, nx; num_boundary_dofs=num_boundary_dofs+1; boundary_dofs(num_boundary_dofs) = ix; end do ! Bottom
    do ix = 1, nx; num_boundary_dofs=num_boundary_dofs+1; boundary_dofs(num_boundary_dofs) = ny*nx+ix; end do ! Top
    do iy = 1, ny; num_boundary_dofs=num_boundary_dofs+1; boundary_dofs(num_boundary_dofs) = num_edges_x+iy; end do ! Left
    do iy = 1, ny; num_boundary_dofs=num_boundary_dofs+1; boundary_dofs(num_boundary_dofs) = num_edges_x+nx*ny+iy; end do ! Right

    call calculate_analytical_dofs_polynomial(U_solution, nx, ny, num_edges_x, dx, dy, 0._PR)
    call calculate_analytical_dofs_polynomial(U_solution_old, nx, ny, num_edges_x, dx, dy,0._PR)

    ! =================================================================
    ! 4. TIME-STEPPING LOOP
    ! =================================================================
    write(*,*) 'Starting time-stepping loop...'
    do time_step = 1, num_time_steps
        current_time = time_step * dt
            ! --- 4a. Reassemble load vector at current time ---
        F_vector = 0.0_PR  ! Reset load vector
        do elem = 1, num_quad_elements
           elem_nodes = quad_connectivity(elem, :)
           do i = 1, 4
              elem_coords(i, :) = nodes(elem_nodes(i), :)
           end do

        ! Edge numbering
            ix = mod(elem-1, nx) + 1
            iy = (elem-1) / nx + 1
            elem_edges(1) = (iy-1)*nx + ix
            elem_edges(2) = num_edges_x + ix*ny + iy
            elem_edges(3) = iy*nx + ix
            elem_edges(4) = num_edges_x + (ix-1)*ny + iy

            ! Call element routine with current_time
            call process_nedelec_element_rhs(elem_nodes, elem_edges, elem_coords, &
                                        Me, Ke, fe, ierr, current_time)
            if (ierr /= 0) stop 'Error in element processing!'

            ! Assemble with sign correction
            s = [1, 1, -1, -1]
            do i = 1, 4
                F_vector(elem_edges(i)) = F_vector(elem_edges(i)) + s(i) * fe(i)
            end do
        end do


        ! --- 4a. Calculate the RHS vector: F_rhs = F_source + (1/dt)*M*U_old ---
        F_rhs = F_vector 
        call sparse_matvec_add(M_sparse, U_solution_old, F_rhs, 1.0_PR / dt)

        ! --- 4b. Apply Dirichlet Boundary Conditions ---
        ! Add debug output
       write(*,*) 'Time step:', time_step, 'Max U_analytical:', maxval(abs(U_analytical))
       write(*,*) 'Time step:', time_step, 'Max U_numerical:', maxval(abs(U_solution))

       !write(*,*) 'F_rhs before BC:', maxval(abs(F_rhs))
       call calculate_analytical_dofs_polynomial(U_analytical, nx, ny, num_edges_x, dx, dy, current_time)
       call apply_dirichlet_bc(A_system_sparse, F_rhs, boundary_dofs, num_boundary_dofs, U_boundary_perfect_conductor)

       !write(*,*) 'F_rhs after BC:', maxval(abs(F_rhs))



       ! print*,maxval(abs(F_rhs))
        call solve_cg(A_system_sparse, F_rhs, U_solution, max_iter, tolerance)



        ! --- 4d. Update for Next Step ---
        U_solution_old = U_solution

    end do
    write(*,*) 'Time-stepping finished.'


    ! =================================================================
    ! 6. ERROR ANALYSIS
    ! =================================================================
   
    call calculate_error(U_solution, U_analytical, num_dofs)

    !print*,maxval(abs(U_solution(boundary_dofs)))



    ! =================================================================
    ! 7. CLEANUP
    ! =================================================================
    deallocate(nodes, quad_connectivity)
    deallocate(U_solution, F_vector, U_analytical, boundary_dofs)
    call sparse_free(K_sparse)
    call sparse_free(M_sparse)
    call triplet_free(K_triplet)
    call triplet_free(M_triplet)
    write(*,*) 'Simulation completed successfully!'

contains

    ! -------------------------------------------------
    ! NEW SUBROUTINES FOR THE STATIC TEST
    ! -------------------------------------------------
subroutine calculate_analytical_dofs_polynomial(U_anal, nx, ny, num_edges_x, dx, dy, current_time)
    implicit none
    real(PR), intent(out) :: U_anal(:)
    integer, intent(in) :: nx, ny, num_edges_x  
    real(PR), intent(in) :: dx, dy, current_time
    integer :: i, ix, iy
    real(PR) :: x_start, x_end, y_start, y_end, x_mid, y_mid
    real(PR) :: Ex_val, Ey_val
    
    U_anal = 0.0_PR
    
    ! Arêtes horizontales (contribution Ex)
    do iy = 1, ny + 1
        do ix = 1, nx
            i = (iy-1)*nx + ix
            
            ! Coordonnées de l'arête horizontale
            x_start = real(ix-1, PR) * dx
            x_end = real(ix, PR) * dx  
            y_mid = real(iy-1, PR) * dy
            x_mid = (x_start + x_end) / 2.0_PR
            
            ! DoF = ∫ Ex dx ≈ Ex(x_mid, y_mid) * dx
            Ex_val = Ax_analytical(x_mid, y_mid, current_time)
            U_anal(i) = Ex_val * dx
        end do
    end do
    
    ! Arêtes verticales (contribution Ey)
    do iy = 1, ny
        do ix = 1, nx + 1
            i = num_edges_x + (ix-1)*ny + iy
            
            ! Coordonnées de l'arête verticale
            x_mid = real(ix-1, PR) * dx
            y_start = real(iy-1, PR) * dy
            y_end = real(iy, PR) * dy
            y_mid = (y_start + y_end) / 2.0_PR
            
            ! DoF = ∫ Ey dy ≈ Ey(x_mid, y_mid) * dy  
            Ey_val = Ay_analytical(x_mid, y_mid, current_time)
            U_anal(i) = Ey_val * dy
        end do
    end do
end subroutine calculate_analytical_dofs_polynomial
subroutine apply_dirichlet_bc(K_mat, F, bnd_dofs, n_bnd, U_anal)
    implicit none
    type(sparse_matrix), intent(inout) :: K_mat
    real(PR), intent(inout) :: F(:)
    integer, intent(in) :: bnd_dofs(:), n_bnd
    real(PR), intent(in) :: U_anal(:)
    
    integer :: i, j, k, row_idx, col_idx
    real(PR) :: val
    logical, allocatable :: is_boundary(:)
    logical, save :: first_call = .true.
    
    ! Storage for original coupling terms (static - persists between calls)
    integer, save :: n_coupling_terms = 0
    integer, allocatable, save :: coupling_rows(:), coupling_cols(:)
    real(PR), allocatable, save :: coupling_values(:)

    ! Create a mask for boundary DOFs
    allocate(is_boundary(K_mat%n))
    is_boundary = .false.
    do i = 1, n_bnd
        is_boundary(bnd_dofs(i)) = .true.
    end do

    ! On first call: store original coupling terms and modify matrix structure
    if (first_call) then
        ! First pass: count coupling terms (interior row, boundary column)
        n_coupling_terms = 0
        do i = 1, K_mat%n
            if (.not. is_boundary(i)) then  ! Interior row
                do k = K_mat%row_ptr(i), K_mat%row_ptr(i+1) - 1
                    col_idx = K_mat%col_indices(k)
                    if (is_boundary(col_idx)) then  ! Boundary column
                        n_coupling_terms = n_coupling_terms + 1
                    end if
                end do
            end if
        end do
        
        ! Allocate storage for coupling terms
        allocate(coupling_rows(n_coupling_terms))
        allocate(coupling_cols(n_coupling_terms))
        allocate(coupling_values(n_coupling_terms))
        
        ! Second pass: store coupling terms
        n_coupling_terms = 0
        do i = 1, K_mat%n
            if (.not. is_boundary(i)) then  ! Interior row
                do k = K_mat%row_ptr(i), K_mat%row_ptr(i+1) - 1
                    col_idx = K_mat%col_indices(k)
                    if (is_boundary(col_idx)) then  ! Boundary column
                        n_coupling_terms = n_coupling_terms + 1
                        coupling_rows(n_coupling_terms) = i
                        coupling_cols(n_coupling_terms) = col_idx
                        coupling_values(n_coupling_terms) = K_mat%values(k)
                    end if
                end do
            end if
        end do
        
        ! Now modify the matrix structure
        ! Zero out rows and columns for boundary DOFs and set diagonal to 1
        do i = 1, n_bnd
            row_idx = bnd_dofs(i)

            ! Zero out the row
            do k = K_mat%row_ptr(row_idx), K_mat%row_ptr(row_idx+1) - 1
                K_mat%values(k) = 0.0_PR
            end do

            ! Zero out the column
            do j = 1, K_mat%n
                do k = K_mat%row_ptr(j), K_mat%row_ptr(j+1) - 1
                    if (K_mat%col_indices(k) == row_idx) then
                        K_mat%values(k) = 0.0_PR
                    end if
                end do
            end do

            ! Set diagonal to 1
            call sparse_set_val(K_mat, row_idx, row_idx, 1.0_PR)
        end do
        
        first_call = .false.
        write(*,*) 'Stored', n_coupling_terms, 'coupling terms for BC application'
    end if

    ! Every call: modify RHS based on stored original coupling values
    do k = 1, n_coupling_terms
        i = coupling_rows(k)
        col_idx = coupling_cols(k)
        val = coupling_values(k)
        F(i) = F(i) - val * U_anal(col_idx)
    end do

    ! Every call: set RHS for boundary DOFs
    do i = 1, n_bnd
        row_idx = bnd_dofs(i)
        F(row_idx) = U_anal(row_idx)
    end do

    deallocate(is_boundary)
end subroutine apply_dirichlet_bc
subroutine calculate_error(U_num, U_anal, n)
        implicit none
        real(PR), intent(in) :: U_num(:), U_anal(:)
        real(PR),allocatable::diff(:)
        integer, intent(in) :: n
        real(PR) :: norm_err, norm_anal

        allocate(diff(n))

            
    
    

        diff = U_num - U_anal
        norm_err = sqrt(dot_product(diff, diff))
        norm_anal = sqrt(dot_product(U_anal, U_anal))

        write(*,*) '========================================='
        write(*,*) 'ERROR ANALYSIS'
        write(*,*) '========================================='
        write(*,'(A, E12.5)') 'L2 Norm of Error ||U_num - U_anal|| = ', norm_err
        write(*,'(A, E12.5)') 'L2 Norm of Analytical ||U_anal||   = ', norm_anal
        if (norm_anal > 1.0e-12_PR) then
            write(*,'(A, E12.5)') 'Relative Error = ', norm_err / norm_anal
        end if
end subroutine calculate_error

end program main
