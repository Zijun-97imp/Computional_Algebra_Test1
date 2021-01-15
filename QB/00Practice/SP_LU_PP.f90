    
#Low-Memory LU factorization

    ! Inputs:
    implicit none
    integer, intent (in) :: n
    integer (4) :: p(1,n)
    real (4) :: A(n,n)
    integer :: k,j,r

    ! Perform the Partial Pivoting Matrix:
    do k = 1,n-1
      !Fine the index of the pivot: r!
      if (A(k,r).NE.0.0) then
        !Perform the rows swaping
        A([k r],:) = A([r k],:)
        p([k r],:) = p([r k],:)
        A(k+1:n,k) = A(k+1:n,k)/A(k,k)
        do j = k+1,n
          A(k+1:n,j) = A(k+1:n,j)-A(k+1:n,k)*A(k,j)
        end do
      else
        print *. "Error: Matrix A is not invertible"
        exit
      end if
    end do