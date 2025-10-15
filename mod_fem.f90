module mod_fem
    use mod_var
    use mod_function
    implicit none
      ! Perméabilité du vide [H/m]

    real(PR), parameter :: mu_0 = 4.0e-7_PR * PI  ! Perméabilité du vide [H/m]
    real(PR), parameter :: sigma_0 = 5.96e7_PR    ! Conductivité du cuivre [S/m]

contains

    ! =================================================================
    ! Subroutine for 4-Node Linear Quadrilateral Elements (Original)
    ! =================================================================


    ! =================================================================
    ! CORRECTED Subroutine for Nédélec First Order Elements on Quadrilaterals
    ! =================================================================
subroutine process_nedelec_element(element_nodes, element_edges, node_coordinates, &
                                       Me, Ke, fe, ierr)
    implicit none
    
    ! Input/Output arguments
    integer, intent(in) :: element_nodes(4)         ! indices of the 4 corner nodes
    integer, intent(in) :: element_edges(4)         ! indices of the 4 edges
    real(PR), intent(in) :: node_coordinates(4,2)   ! (x,y) coordinates of nodes
    real(PR), intent(out) :: Ke(4,4), Me(4,4)       ! element stiffness and mass matrices
    real(PR), intent(out) :: fe(4)                  ! element load vector
    integer, intent(out) :: ierr                    ! error flag
    
    ! Local variables
    real(PR) :: x(4), y(4)
    real(PR) :: gauss_points(2), weights(2)
    real(PR) :: xi, eta, w_xi, w_eta
    real(PR) :: N_edge(4,2)                         ! Nédélec shape functions (vectorial) in physical space
    real(PR) :: curl_N_edge(4)                      ! Curl of Nédélec shape functions in physical space
    real(PR) :: dx_dxi, dx_deta, dy_dxi, dy_deta
    real(PR) :: det_J, inv_det_J
    real(PR) :: J_inv_11, J_inv_12, J_inv_21, J_inv_22
    real(PR) :: x_pt, y_pt
    real(PR) :: sigma, mu_r, weight
    integer :: i, j, m, n
    
    ! Shape function derivatives for coordinate transformation
    real(PR) :: dN_dxi(4), dN_deta(4)
    real(PR) :: N_nodal(4)
    
    ! Reference edge functions and their curls
    real(PR) :: ref_N_edge(4,2)
    real(PR) :: ref_curl_N_edge(4)
    
    ! Initialize
    ierr = 0
    
    ! Extract coordinates
    do i = 1, 4
        x(i) = node_coordinates(i,1)
        y(i) = node_coordinates(i,2)
    end do
    
    ! Set up integration (2x2 Gauss points)
    gauss_points(1) = -0.5773502692_PR  ! -1/sqrt(3)
    gauss_points(2) =  0.5773502692_PR  ! +1/sqrt(3)
    weights(1) = 1.0_PR
    weights(2) = 1.0_PR
    
    ! Initialize element matrices
    Ke = 0.0_PR
    fe = 0.0_PR
    Me = 0.0_PR
    
    ! Curls of the reference Nédélec functions (constant for the lowest order)

    ref_curl_N_edge(1) =  0.25_PR
    ref_curl_N_edge(2) =  0.25_PR
    ref_curl_N_edge(3) =  0.25_PR
    ref_curl_N_edge(4) =  0.25_PR

    ! Loop over integration points
    do i = 1, 2
        do j = 1, 2
            xi = gauss_points(i)
            eta = gauss_points(j)
            w_xi = weights(i)
            w_eta = weights(j)

            ! Define reference Nédélec functions at the current integration point (xi, eta)
            ! These are vector functions, not constants.
            ref_N_edge(1,1) = 0.25_PR * (1.0_PR - eta)
            ref_N_edge(1,2) = 0.0_PR

            ref_N_edge(2,1) = 0.0_PR
            ref_N_edge(2,2) = 0.25_PR * (1.0_PR + xi)

            ref_N_edge(3,1) = -0.25_PR * (1.0_PR + eta)
            ref_N_edge(3,2) = 0.0_PR

            ref_N_edge(4,1) = 0.0_PR
            ref_N_edge(4,2) = -0.25_PR * (1.0_PR - xi)
            
            ! Compute nodal shape functions for coordinate transformation
            N_nodal(1) = (1.0_PR - xi) * (1.0_PR - eta) / 4.0_PR
            N_nodal(2) = (1.0_PR + xi) * (1.0_PR - eta) / 4.0_PR
            N_nodal(3) = (1.0_PR + xi) * (1.0_PR + eta) / 4.0_PR
            N_nodal(4) = (1.0_PR - xi) * (1.0_PR + eta) / 4.0_PR
            
            ! Derivatives of nodal shape functions
            dN_dxi(1) = -(1.0_PR - eta) / 4.0_PR
            dN_dxi(2) = +(1.0_PR - eta) / 4.0_PR
            dN_dxi(3) = +(1.0_PR + eta) / 4.0_PR
            dN_dxi(4) = -(1.0_PR + eta) / 4.0_PR
            
            dN_deta(1) = -(1.0_PR - xi) / 4.0_PR
            dN_deta(2) = -(1.0_PR + xi) / 4.0_PR
            dN_deta(3) = +(1.0_PR + xi) / 4.0_PR
            dN_deta(4) = +(1.0_PR - xi) / 4.0_PR
            
            ! Compute Jacobian matrix components
            dx_dxi = dot_product(dN_dxi, x)
            dx_deta = dot_product(dN_deta, x)
            dy_dxi = dot_product(dN_dxi, y)
            dy_deta = dot_product(dN_deta, y)
            
            ! Compute Jacobian determinant
            det_J = dx_dxi * dy_deta - dx_deta * dy_dxi
           ! print *, "det_J = ", det_J
            
            if (det_J <= 0.0_PR) then
                ierr = 1
                return
            end if
            
            ! Compute inverse Jacobian
            inv_det_J = 1.0_PR / det_J
            J_inv_11 =  dy_deta * inv_det_J
            J_inv_12 = -dy_dxi * inv_det_J
            J_inv_21 = -dx_deta * inv_det_J
            J_inv_22 =  dx_dxi * inv_det_J
            
            ! Transform Nédélec functions to physical space using Piola transform
            ! N_phys = (J^-T) * N_ref
            do n = 1, 4
                N_edge(n,1) = J_inv_11 * ref_N_edge(n,1) + J_inv_21 * ref_N_edge(n,2)
                N_edge(n,2) = J_inv_12 * ref_N_edge(n,1) + J_inv_22 * ref_N_edge(n,2)
            end do
            
            ! Transform curl to physical space
            ! curl_phys = curl_ref / det_J
            do n = 1, 4
                curl_N_edge(n) = ref_curl_N_edge(n) / det_J
            end do
            
            ! Compute physical coordinates at integration point
            x_pt = dot_product(N_nodal, x)
            y_pt = dot_product(N_nodal, y)
            
            ! Evaluate material properties
            sigma = sigma_function(x_pt, y_pt)     ! Conductivity
            mu_r = mu_r_function(x_pt, y_pt)       ! Relative permeability
            
            ! Compute integration weight
            weight = w_xi * w_eta * det_J
            
            ! Assemble element stiffness matrix (curl-curl term)
            ! Ke = ∫ (1/μ) (curl w_a) · (curl w_b) dΩ
            do m = 1, 4
                do n = 1, 4
                    Ke(m,n) = Ke(m,n) + weight * (1.0_PR/mu_r) * &
                              curl_N_edge(m) * curl_N_edge(n)
                end do
            end do
            
            ! Assemble element mass matrix (conductivity term)
            ! Me = ∫ σ (w_a · w_b) dΩ
            do m = 1, 4
                do n = 1, 4
                    Me(m,n) = Me(m,n) + weight * sigma * &
                              (N_edge(m,1)*N_edge(n,1) + N_edge(m,2)*N_edge(n,2))
                end do
            end do
            
            ! Assemble element load vector
            ! fe = ∫ (J_s · w_a) dΩ
            !do m = 1, 4
                !fe(m) = fe(m) + weight * (Js_x_function(x_pt, y_pt,) * N_edge(m,1) + &
                                          !Js_y_function(x_pt, y_pt) * N_edge(m,2))
            !end do
            
        end do
    end do
    
end subroutine process_nedelec_element

subroutine process_nedelec_element_rhs(element_nodes, element_edges, node_coordinates, &
                                       Me, Ke, fe, ierr,time)
    implicit none
    
    ! Input/Output arguments
    integer, intent(in) :: element_nodes(4)         ! indices of the 4 corner nodes
    integer, intent(in) :: element_edges(4)         ! indices of the 4 edges
    real(PR), intent(in) :: node_coordinates(4,2) 
    real(PR), intent(in) :: time   ! (x,y) coordinates of nodes
    real(PR), intent(out) :: Ke(4,4), Me(4,4)       ! element stiffness and mass matrices
    real(PR), intent(out) :: fe(4)                  ! element load vector
    integer, intent(out) :: ierr                    ! error flag
    
    ! Local variables
    real(PR) :: x(4), y(4)
    real(PR) :: gauss_points(2), weights(2)
    real(PR) :: xi, eta, w_xi, w_eta
    real(PR) :: N_edge(4,2)                         ! Nédélec shape functions (vectorial) in physical space
    real(PR) :: curl_N_edge(4)                      ! Curl of Nédélec shape functions in physical space
    real(PR) :: dx_dxi, dx_deta, dy_dxi, dy_deta
    real(PR) :: det_J, inv_det_J
    real(PR) :: J_inv_11, J_inv_12, J_inv_21, J_inv_22
    real(PR) :: x_pt, y_pt
    real(PR) :: sigma, mu_r, weight
    integer :: i, j, m, n
    
    ! Shape function derivatives for coordinate transformation
    real(PR) :: dN_dxi(4), dN_deta(4)
    real(PR) :: N_nodal(4)
    
    ! Reference edge functions and their curls
    real(PR) :: ref_N_edge(4,2)
    real(PR) :: ref_curl_N_edge(4)
    
    ! Initialize
    ierr = 0
    
    ! Extract coordinates
    do i = 1, 4
        x(i) = node_coordinates(i,1)
        y(i) = node_coordinates(i,2)
    end do
    
    ! Set up integration (2x2 Gauss points)
    gauss_points(1) = -0.5773502692_PR  ! -1/sqrt(3)
    gauss_points(2) =  0.5773502692_PR  ! +1/sqrt(3)
    weights(1) = 1.0_PR
    weights(2) = 1.0_PR
    
    ! Initialize element matrices
    Ke = 0.0_PR
    fe = 0.0_PR
    Me = 0.0_PR
    
    ! Curls of the reference Nédélec functions (constant for the lowest order)

    ref_curl_N_edge(1) =  0.25_PR
    ref_curl_N_edge(2) =  0.25_PR
    ref_curl_N_edge(3) =  0.25_PR
    ref_curl_N_edge(4) =  0.25_PR

    ! Loop over integration points
    do i = 1, 2
        do j = 1, 2
            xi = gauss_points(i)
            eta = gauss_points(j)
            w_xi = weights(i)
            w_eta = weights(j)

            ! Define reference Nédélec functions at the current integration point (xi, eta)
            ! These are vector functions, not constants.
            ref_N_edge(1,1) = 0.25_PR * (1.0_PR - eta)
            ref_N_edge(1,2) = 0.0_PR

            ref_N_edge(2,1) = 0.0_PR
            ref_N_edge(2,2) = 0.25_PR * (1.0_PR + xi)

            ref_N_edge(3,1) = -0.25_PR * (1.0_PR + eta)
            ref_N_edge(3,2) = 0.0_PR

            ref_N_edge(4,1) = 0.0_PR
            ref_N_edge(4,2) = -0.25_PR * (1.0_PR - xi)
            
            ! Compute nodal shape functions for coordinate transformation
            N_nodal(1) = (1.0_PR - xi) * (1.0_PR - eta) / 4.0_PR
            N_nodal(2) = (1.0_PR + xi) * (1.0_PR - eta) / 4.0_PR
            N_nodal(3) = (1.0_PR + xi) * (1.0_PR + eta) / 4.0_PR
            N_nodal(4) = (1.0_PR - xi) * (1.0_PR + eta) / 4.0_PR
            
            ! Derivatives of nodal shape functions
            dN_dxi(1) = -(1.0_PR - eta) / 4.0_PR
            dN_dxi(2) = +(1.0_PR - eta) / 4.0_PR
            dN_dxi(3) = +(1.0_PR + eta) / 4.0_PR
            dN_dxi(4) = -(1.0_PR + eta) / 4.0_PR
            
            dN_deta(1) = -(1.0_PR - xi) / 4.0_PR
            dN_deta(2) = -(1.0_PR + xi) / 4.0_PR
            dN_deta(3) = +(1.0_PR + xi) / 4.0_PR
            dN_deta(4) = +(1.0_PR - xi) / 4.0_PR
            
            ! Compute Jacobian matrix components
            dx_dxi = dot_product(dN_dxi, x)
            dx_deta = dot_product(dN_deta, x)
            dy_dxi = dot_product(dN_dxi, y)
            dy_deta = dot_product(dN_deta, y)
            
            ! Compute Jacobian determinant
            det_J = dx_dxi * dy_deta - dx_deta * dy_dxi
           ! print *, "det_J = ", det_J
            
            if (det_J <= 0.0_PR) then
                ierr = 1
                return
            end if
            
            ! Compute inverse Jacobian
            inv_det_J = 1.0_PR / det_J
            J_inv_11 =  dy_deta * inv_det_J
            J_inv_12 = -dy_dxi * inv_det_J
            J_inv_21 = -dx_deta * inv_det_J
            J_inv_22 =  dx_dxi * inv_det_J
            
            ! Transform Nédélec functions to physical space using Piola transform
            ! N_phys = (J^-T) * N_ref
            do n = 1, 4
                N_edge(n,1) = J_inv_11 * ref_N_edge(n,1) + J_inv_21 * ref_N_edge(n,2)
                N_edge(n,2) = J_inv_12 * ref_N_edge(n,1) + J_inv_22 * ref_N_edge(n,2)
            end do
            
            ! Transform curl to physical space
            ! curl_phys = curl_ref / det_J
            do n = 1, 4
                curl_N_edge(n) = ref_curl_N_edge(n) / det_J
            end do
            
            ! Compute physical coordinates at integration point
            x_pt = dot_product(N_nodal, x)
            y_pt = dot_product(N_nodal, y)
            
            ! Evaluate material properties
            sigma = sigma_function(x_pt, y_pt)     ! Conductivity
            mu_r = mu_r_function(x_pt, y_pt)       ! Relative permeability
            
            ! Compute integration weight
            weight = w_xi * w_eta * det_J

            
            ! Assemble element load vector
            ! fe = ∫ (J_s · w_a) dΩ
            do m = 1, 4
                fe(m) = fe(m) + weight * (Js_x_function(x_pt, y_pt,time) * N_edge(m,1) + &
                                          Js_y_function(x_pt, y_pt,time) * N_edge(m,2))
            end do
            
        end do
    end do
    
end subroutine process_nedelec_element_rhs
    
    ! =================================================================
    ! CORRECTED function to find boundary edges
    ! =================================================================
    subroutine find_boundary_edges(nx, ny, num_edges_x, num_edges_y, boundary_dofs, num_boundary_dofs)
        implicit none
        integer, intent(in) :: nx, ny, num_edges_x, num_edges_y
        integer, intent(out) :: boundary_dofs(:)
        integer, intent(out) :: num_boundary_dofs
        integer :: ix, iy
        
        num_boundary_dofs = 0
        
        ! Arêtes horizontales du bord inférieur (y=0)
        do ix = 1, nx
            num_boundary_dofs = num_boundary_dofs + 1
            boundary_dofs(num_boundary_dofs) = ix
        end do
        
        ! Arêtes horizontales du bord supérieur (y=Ly)
        do ix = 1, nx
            num_boundary_dofs = num_boundary_dofs + 1
            boundary_dofs(num_boundary_dofs) = ny * nx + ix
        end do
        
        ! Arêtes verticales du bord gauche (x=0)
        do iy = 1, ny
            num_boundary_dofs = num_boundary_dofs + 1
            boundary_dofs(num_boundary_dofs) = num_edges_x + iy
        end do
        
        ! Arêtes verticales du bord droit (x=Lx)
        do iy = 1, ny
            num_boundary_dofs = num_boundary_dofs + 1
            boundary_dofs(num_boundary_dofs) = num_edges_x + nx * ny + iy
        end do
        
    end subroutine find_boundary_edges
    
    ! =================================================================
    ! Auxiliary functions for magnetodynamics
    ! =================================================================
    
    function sigma_function(x, y) result(sigma)
        implicit none
        real(PR), intent(in) :: x, y
        real(PR) :: sigma
        ! Conductivity function - can be spatially varying
        sigma = 1.0_PR  ! Default: uniform conductivity
        
        ! Example: conductivity varies with position
        ! sigma = 1.0_PR + 0.5_PR * sin(2.0_PR * 3.14159_PR * x) * sin(2.0_PR * 3.14159_PR * y)
    end function sigma_function
    
    function mu_r_function(x, y) result(mu_r)
        implicit none
        real(PR), intent(in) :: x, y
        real(PR) :: mu_r
        ! Relative permeability function
        mu_r = 1.0_PR  ! Default: vacuum permeability
        
        ! Example: magnetic material in center
        ! if ((x-0.5_PR)**2 + (y-0.5_PR)**2 < 0.1_PR) then
        !     mu_r = 1000.0_PR  ! High permeability material
        ! else
        !     mu_r = 1.0_PR
        ! end if
    end function mu_r_function
    
    function Js_x_function(x, y, time) result(Js_x)
        implicit none
        real(PR), intent(in) :: x, y, time
        real(PR) :: Js_x
        real(PR), parameter :: PI = 3.14159265358979323846_PR
        real(PR) :: mu_r, sigma, factor

        ! J = (2*pi^2/mu - sigma) * A
        mu_r = mu_r_function(x, y)
        sigma = sigma_function(x, y)
        factor = (2.0_PR * PI**2 / mu_r - sigma)
        Js_x = factor * Ax_analytical(x, y, time)
    end function Js_x_function

    function Js_y_function(x, y, time) result(Js_y)
        implicit none
        real(PR), intent(in) :: x, y, time
        real(PR) :: Js_y
        real(PR), parameter :: PI = 3.14159265358979323846_PR
        real(PR) :: mu_r, sigma, factor

        ! J = (2*pi^2/mu - sigma) * A
        mu_r = mu_r_function(x, y)
        sigma = sigma_function(x, y)
        factor = (2.0_PR * PI**2 /mu_r - sigma)
        Js_y = factor * Ay_analytical(x, y, time)
    end function Js_y_function

    ! =================================================================
    ! Utility function to set initial conditions for magnetodynamics
    ! =================================================================
    subroutine set_initial_field(U_n, nx, ny, num_edges_x, Lx, Ly, dx, dy)
        implicit none
        real(PR), intent(out) :: U_n(:)
        integer, intent(in) :: nx, ny, num_edges_x
        real(PR), intent(in) :: Lx, Ly, dx, dy
        
        integer :: ix, iy, i
        real(PR) :: x, y
        real(PR), parameter :: PI = 3.14159265358979323846_PR
        
        U_n = 0.0_PR
        
        ! Set initial field - example: decaying vortex
        ! Horizontal edges (Ex component)
        do iy = 1, ny+1
            do ix = 1, nx
                i = (iy-1)*nx + ix
                x = (ix - 0.5_PR) * dx
                y = (iy - 1.0_PR) * dy
                ! Vortex field - Ex component
                U_n(i) = -(y - Ly/2.0_PR) * exp(-10.0_PR * ((x-Lx/2.0_PR)**2 + (y-Ly/2.0_PR)**2))
            end do
        end do
        
        ! Vertical edges (Ey component)  
        do iy = 1, ny
            do ix = 1, nx+1
                i = num_edges_x + (ix-1)*ny + iy
                x = (ix - 1.0_PR) * dx
                y = (iy - 0.5_PR) * dy
                ! Vortex field - Ey component
                U_n(i) = (x - Lx/2.0_PR) * exp(-10.0_PR * ((x-Lx/2.0_PR)**2 + (y-Ly/2.0_PR)**2))
            end do
        end do
        
    end subroutine set_initial_field

  

    ! Basic linear solver (kept for compatibility)
   

    ! =================================================================
    ! Analytical Solution for Verification (TM11 mode in a cavity)
    ! =================================================================

   


end module mod_fem
