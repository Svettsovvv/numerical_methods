program task_one
    implicit none
    integer, parameter :: n = 20
    real :: A(n, n), B(n), X(n)
    real :: P(3) = [1.0000000, 0.0156250, 0.0000000]
    integer :: i, j, k
    real :: temp

    open(unit=10, file='RGEmtr094.dat', status='old')
    do i = 1, n
        read(10, *) (A(i, j), j = 1, n)
    end do
    close(10)

    open(unit=11, file='RGErhs094.dat', status='old')
    do i = 1, n
        read(11, *) B(i)
    end do
    close(11)

    do k = 1, 3
        A(1, 1) = A(1, 1) + P(k)
        B(1) = B(1) + 8 * P(k)
    end do

    call gaussian_method(n, A, B, X)

    print *, 'Solve system:'
    do i = 1, n
        print *, X(i)
    end do

contains

    subroutine gaussian_method(n, A, B, X)
        integer, intent(in) :: n
        real, intent(inout) :: A(n, n), B(n)
        real, intent(out) :: X(n)
        integer :: i, j, k
        real :: factor

        do k = 1, n-1
            do i = k+1, n
                factor = A(i, k) / A(k, k)
                do j = k+1, n
                    A(i, j) = A(i, j) - factor * A(k, j)
                end do
                B(i) = B(i) - factor * B(k)
            end do
        end do

        X(n) = B(n) / A(n, n)
        do i = n-1, 1, -1
            X(i) = B(i)
            do j = i+1, n
                X(i) = X(i) - A(i, j) * X(j)
            end do
            X(i) = X(i) / A(i, i)
        end do
    end subroutine gaussian_method

end program task_one
