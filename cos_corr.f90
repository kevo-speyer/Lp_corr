subroutine r_next_corr(R_next, R_mean, R_var, n_mon_1, n_dim, corr_out)
implicit none
integer :: i, j, i_dim
integer, intent(in) :: n_mon_1, n_dim
real(kind=8), dimension(n_mon_1,n_dim), intent(in) :: R_next, R_mean
real(kind=8), dimension(n_mon_1), intent(in) :: R_var
real(kind=8), dimension(0:n_mon_1-1), intent(out) :: corr_out
real(kind=8), dimension(0:n_mon_1-1) :: n_event
real(kind=8) :: corr_dummy

corr_out(:) = 0.
n_event(:) = 0.

do i = 1, n_mon_1 !loop  over monomes 
    do j = i, n_mon_1 ! loop over monomers
      
        n_event(j-i) = n_event(j-i) + 1
        corr_dummy = 0.
       
        do i_dim = 1, 3
            corr_dummy = corr_dummy + (R_next(i,i_dim) - R_mean(i,i_dim)) * (R_next(j,i_dim) - R_mean(j,i_dim))          
        end do
     
        corr_out(j-i) = corr_out(j-i) + corr_dummy/sqrt( R_var(i) * R_var(j))

    end do
end do

!Normalization
corr_out(:) = corr_out(:) / n_event(:)

!corr_out(:) = (corr_out(:) - meanR_sq) / (mean_r2 - meanR_sq)

end subroutine

subroutine cos_corr(R_next, n_mon_1, n_dim, corr_out)
implicit none
integer :: i, j, i_dim
integer, intent(in) :: n_mon_1, n_dim
real(kind=8), dimension(n_dim,n_mon_1), intent(in) :: R_next
real(kind=8), dimension(0:n_mon_1-1), intent(out) :: corr_out
real(kind=8), dimension(0:n_mon_1-1) :: n_event
real(kind=8) :: dot_prod, norm_r2, norm_r1

corr_out(:) = 0.
n_event(:) = 0.

do i = 1, n_mon_1 !loop  over monomes 
    do j = i, n_mon_1 ! loop over monomers
      
        n_event(j-i) = n_event(j-i) + 1
        dot_prod = 0.
        norm_r1 = 0.
        norm_r2 = 0.
       
        do i_dim = 1, 3
            dot_prod = dot_prod + R_next(i_dim,i) * R_next(i_dim,j)       
            norm_r1 = norm_r1 +  R_next(i_dim,i) * R_next(i_dim,i)   
            norm_r2 = norm_r2 +  R_next(i_dim,j) * R_next(i_dim,j)
        end do
     
        corr_out(j-i) = corr_out(j-i) + dot_prod / sqrt( norm_r1 * norm_r2 )

    end do
end do

!Normalization
corr_out(:) = corr_out(:) / n_event(:)

!corr_out(:) = (corr_out(:) - meanR_sq) / (mean_r2 - meanR_sq)

end subroutine 
