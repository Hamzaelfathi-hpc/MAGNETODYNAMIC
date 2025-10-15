module mod_function
    use mod_var
    implicit none

    type :: sparse_matrix
        integer :: n                      ! Matrix dimension (n x n)
        integer :: nnz                    ! Number of non-zero elements
        real(PR), allocatable :: values(:)      ! Non-zero values
        integer, allocatable :: col_indices(:)! Column indices
        integer, allocatable :: row_ptr(:)      ! Row pointers
    end type sparse_matrix

    ! Triplet format for assembly
    type :: triplet_matrix
        integer :: nrows, ncols, nnz, capacity
        integer, allocatable :: row(:), col(:)
        real(PR), allocatable :: val(:)
    end type triplet_matrix

contains
subroutine sparse_set_val(K_mat, row, col, val)
    type(sparse_matrix), intent(inout) :: K_mat
    integer, intent(in) :: row, col
    real(PR), intent(in) :: val
    
    integer :: k, pos
    integer, allocatable :: temp_cols(:)
    real(PR), allocatable :: temp_vals(:)
    logical :: found
    
    ! Check if entry exists
    found = .false.
    do k = K_mat%row_ptr(row), K_mat%row_ptr(row+1)-1
        if (K_mat%col_indices(k) == col) then
            K_mat%values(k) = val
            found = .true.
            exit
        end if
    end do
    
    ! If not found, insert new entry
    if (.not. found) then
        ! Expand storage if needed
        if (K_mat%nnz + 1 > size(K_mat%values)) then
            allocate(temp_cols(size(K_mat%col_indices)))
            allocate(temp_vals(size(K_mat%values)))
            temp_cols = K_mat%col_indices
            temp_vals = K_mat%values
            deallocate(K_mat%col_indices, K_mat%values)
            
            allocate(K_mat%col_indices(size(temp_cols) + 100))
            allocate(K_mat%values(size(temp_vals) + 100))
            K_mat%col_indices(1:size(temp_cols)) = temp_cols
            K_mat%values(1:size(temp_vals)) = temp_vals
            deallocate(temp_cols, temp_vals)
        end if
        
        ! Make space for new entry
        do pos = K_mat%nnz, K_mat%row_ptr(row+1), -1
            K_mat%col_indices(pos+1) = K_mat%col_indices(pos)
            K_mat%values(pos+1) = K_mat%values(pos)
        end do
        
        ! Insert new entry
        K_mat%col_indices(K_mat%row_ptr(row+1)) = col
        K_mat%values(K_mat%row_ptr(row+1)) = val
        K_mat%nnz = K_mat%nnz + 1
        
        ! Update row pointers
        do k = row+1, K_mat%n+1
            K_mat%row_ptr(k) = K_mat%row_ptr(k) + 1
        end do
    end if
end subroutine sparse_set_val

    ! Initialize triplet matrix
    subroutine triplet_init(T, nrows, estimated_nnz)
        type(triplet_matrix), intent(out) :: T
        integer, intent(in) :: nrows, estimated_nnz
        
        T%nrows = nrows
        T%ncols = nrows
        T%nnz = 0
        T%capacity = estimated_nnz
        
        allocate(T%row(estimated_nnz))
        allocate(T%col(estimated_nnz))
        allocate(T%val(estimated_nnz))
    end subroutine

    ! Add entry to triplet matrix by appending. Duplicates will be handled later.
    subroutine triplet_add_entry(T, i, j, value)
        type(triplet_matrix), intent(inout) :: T
        integer, intent(in) :: i, j
        real(PR), intent(in) :: value
        
        integer, allocatable :: temp_row(:), temp_col(:)
        real(PR), allocatable :: temp_val(:)
        integer :: new_capacity
        
        if (abs(value) < 1.0e-14_PR) return  ! Skip near-zero values
        
        ! Expand storage if needed
        if (T%nnz >= T%capacity) then
            new_capacity = T%capacity * 2
            if (new_capacity == 0) new_capacity = 16 ! Initial capacity
            
            allocate(temp_row(new_capacity))
            allocate(temp_col(new_capacity))
            allocate(temp_val(new_capacity))
            
            if (T%nnz > 0) then
                temp_row(1:T%nnz) = T%row(1:T%nnz)
                temp_col(1:T%nnz) = T%col(1:T%nnz)
                temp_val(1:T%nnz) = T%val(1:T%nnz)
            end if
            
            call move_alloc(temp_row, T%row)
            call move_alloc(temp_col, T%col)
            call move_alloc(temp_val, T%val)
            
            T%capacity = new_capacity
        end if
        
        ! Append the new entry. No searching for duplicates.
        T%nnz = T%nnz + 1
        T%row(T%nnz) = i
        T%col(T%nnz) = j
        T%val(T%nnz) = value
    end subroutine

    ! Convert triplet list (with potential duplicates) to CSR format.
    subroutine triplet_to_csr(T, S)
        type(triplet_matrix), intent(in) :: T
        type(sparse_matrix), intent(out) :: S
        
        integer :: i, j, k, p, p_prev
        integer, allocatable :: perm(:)
        integer :: final_nnz
        real(PR), allocatable :: temp_v(:)
        integer, allocatable :: temp_c(:)

        S%n = T%nrows
        S%nnz = 0

        if (T%nnz == 0) then
            allocate(S%row_ptr(S%n + 1))
            S%row_ptr = 1
            return
        end if
        
        allocate(perm(T%nnz))
        do k = 1, T%nnz
            perm(k) = k
        end do
        
        do i = 2, T%nnz
            j = i
            do while (j > 1)
                p_prev = perm(j-1)
                p = perm(j)
                if (T%row(p) < T%row(p_prev) .or. &
                   (T%row(p) == T%row(p_prev) .and. T%col(p) < T%col(p_prev))) then
                    k = perm(j)
                    perm(j) = perm(j-1)
                    perm(j-1) = k
                    j = j - 1
                else
                    exit
                end if
            end do
        end do
        
        allocate(S%values(T%nnz), S%col_indices(T%nnz))
        allocate(S%row_ptr(S%n + 1))
        
        p = perm(1)
        final_nnz = 1
        S%values(final_nnz) = T%val(p)
        S%col_indices(final_nnz) = T%col(p)
        
        do i = 1, T%row(p)
           S%row_ptr(i) = 1
        end do

        do k = 2, T%nnz
            p = perm(k)
            p_prev = perm(k-1)
            
            if (T%row(p) == T%row(p_prev) .and. T%col(p) == T%col(p_prev)) then
                S%values(final_nnz) = S%values(final_nnz) + T%val(p)
            else
                if (T%row(p) > T%row(p_prev)) then
                   do i = T%row(p_prev) + 1, T%row(p)
                      S%row_ptr(i) = final_nnz + 1
                   end do
                end if

                final_nnz = final_nnz + 1
                S%values(final_nnz) = T%val(p)
                S%col_indices(final_nnz) = T%col(p)
            end if
        end do
        
        p = perm(T%nnz)
        do i = T%row(p) + 1, S%n + 1
            S%row_ptr(i) = final_nnz + 1
        end do
        
        S%nnz = final_nnz
        
        if (S%nnz > 0) then
            allocate(temp_v(S%nnz), source=S%values(1:S%nnz))
            call move_alloc(temp_v, S%values)
            allocate(temp_c(S%nnz), source=S%col_indices(1:S%nnz))
            call move_alloc(temp_c, S%col_indices)
        else
            if (allocated(S%values)) deallocate(S%values)
            if (allocated(S%col_indices)) deallocate(S%col_indices)
        end if
    end subroutine
    
    ! Free triplet matrix
    subroutine triplet_free(T)
        type(triplet_matrix), intent(inout) :: T
        
        if (allocated(T%row)) deallocate(T%row)
        if (allocated(T%col)) deallocate(T%col)
        if (allocated(T%val)) deallocate(T%val)
        T%nnz = 0
        T%capacity = 0
    end subroutine

    ! Free sparse matrix
    subroutine sparse_free(S)
        type(sparse_matrix), intent(inout) :: S
        
        if (allocated(S%row_ptr)) deallocate(S%row_ptr)
        if (allocated(S%col_indices)) deallocate(S%col_indices)
        if (allocated(S%values)) deallocate(S%values)
        S%nnz = 0
    end subroutine
    
    ! Standard sparse matrix-vector multiplication
    subroutine sparse_matvec(A, x, y)
        type(sparse_matrix), intent(in) :: A
        real(PR), intent(in) :: x(:)
        real(PR), intent(out) :: y(:)
        
        integer :: i, j, start_idx, end_idx
        
        y = 0.0_PR
        
        do i = 1, A%n
            start_idx = A%row_ptr(i)
            end_idx = A%row_ptr(i+1) - 1
            
            do j = start_idx, end_idx
                y(i) = y(i) + A%values(j) * x(A%col_indices(j))
            end do
        end do
    end subroutine



    ! Material property functions
    function epsilon_function(x, y) result(eps)
        real(PR), intent(in) :: x, y
        real(PR) :: eps
        eps = 1.0_PR
    end function epsilon_function

    function rho_function(x, y) result(rho)
        real(PR), intent(in) :: x, y
        real(PR) :: rho
        rho = 1.0_PR
    end function rho_function

    subroutine solve_cg(K_sparse, F_vector, U_solution, max_iter, tolerance)
    ! Solves the sparse linear system K*U = F using the Conjugate Gradient method.
    ! Assumes K is symmetric and positive-definite.
    implicit none

    ! Arguments
    type(sparse_matrix), intent(in) :: K_sparse
    real(PR), intent(in)            :: F_vector(:)
    real(PR), intent(inout)         :: U_solution(:)  ! Initial guess on input, solution on output
    integer, intent(in)             :: max_iter
    real(PR), intent(in)             :: tolerance

    ! Local variables
    integer :: n, iter
    real(PR), allocatable :: r(:), p(:), Ap(:)
    real(PR) :: rs_old, rs_new, alpha

    n = size(F_vector)
    allocate(r(n), p(n), Ap(n))

    ! 1. Initialize
    ! r = F - K*U  (initial residual)
    call sparse_matvec(K_sparse, U_solution, Ap)
    r = F_vector - Ap
    
    ! p = r  (initial search direction)
    p = r
    
    ! rs_old = r' * r
    rs_old = dot_product(r, r)

    if (sqrt(rs_old) < tolerance) then
        write(*,*) 'CG: Initial guess is already a solution.'
        deallocate(r, p, Ap)
        return
    end if

    ! 2. Main Iteration Loop
    do iter = 1, max_iter
        ! Ap = K * p
        call sparse_matvec(K_sparse, p, Ap)
        
        ! alpha = (r' * r) / (p' * Ap)
        alpha = rs_old / dot_product(p, Ap)
        
        ! U = U + alpha * p (update solution)
        U_solution = U_solution + alpha * p
        
        ! r = r - alpha * Ap (update residual)
        r = r - alpha * Ap
        
        ! rs_new = r' * r
        rs_new = dot_product(r, r)
        
        ! Check for convergence
        if (sqrt(rs_new) < tolerance) then
            !write(*,*) 'CG: Converged in', iter, 'iterations.'
            deallocate(r, p, Ap)
            return
        end if
        
        ! p = r + (rs_new / rs_old) * p (update search direction)
        p = r + (rs_new / rs_old) * p
        
        ! rs_old = rs_new
        rs_old = rs_new
    end do

    write(*,*) 'CG: Reached maximum iterations without converging.'
    deallocate(r, p, Ap)

end subroutine solve_cg

    subroutine add_scaled_csr_to_triplet(trip, csr, scale)
        ! Helper to add all entries from a scaled CSR matrix to a triplet matrix
        implicit none
        type(triplet_matrix), intent(inout) :: trip
        type(sparse_matrix), intent(in) :: csr
        real(PR), intent(in) :: scale
        integer :: i, k, row

        row = 1
        do k = 1, csr%nnz
            do while (k >= csr%row_ptr(row+1))
                row = row + 1
            end do
            call triplet_add_entry(trip, row, csr%col_indices(k), csr%values(k) * scale)
        end do
    end subroutine add_scaled_csr_to_triplet

    subroutine sparse_matvec_add(A, x, y, alpha)
        ! Computes y = y + alpha * A*x
        implicit none
        type(sparse_matrix), intent(in) :: A
        real(PR), intent(in) :: x(:)
        real(PR), intent(inout) :: y(:)
        real(PR), intent(in) :: alpha
        integer :: i, k
        do i = 1, A%n
            do k = A%row_ptr(i), A%row_ptr(i+1) - 1
                y(i) = y(i) + alpha * A%values(k) * x(A%col_indices(k))
            end do
        end do
    end subroutine sparse_matvec_add

function Ax_analytical(x, y, t) result(Ax)
    implicit none
    real(PR), intent(in) :: x, y, t
    real(PR) :: Ax
    real(PR), parameter :: PI = 3.14159265358979323846_PR

    ! Manufactured solution for Ax that is zero at y=0 and y=1
    Ax = sin(PI * y) * cos(PI * x) * exp(-t)
end function Ax_analytical

function Ay_analytical(x, y, t) result(Ay)
    implicit none
    real(PR), intent(in) :: x, y, t
    real(PR) :: Ay
    real(PR), parameter :: PI = 3.14159265358979323846_PR

    ! Manufactured solution for Ay that is zero at x=0 and x=1
    Ay = -cos(PI * y) * sin(PI * x) * exp(-t)
end function Ay_analytical





end module mod_function
