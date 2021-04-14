!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++
! copied from selfe2netcdf interpolatez subroutine 
! Please keep in sync!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!Generate wrapper by calling
!f2py -m interpz -h interpz.pyf interpz.f90
!Generate library
!f2py  -c interpz.pyf  interpz.f90

subroutine interpz1d(varin, zin, zout, np, nzin, nzout, kz, null_value, varout)


    integer                     , intent(in)        :: np, nzin, nzout, kz
    real(8)                     , intent(in)        :: null_value
    real(8), dimension(np,nzin) , intent(in)        :: zin , varin
    real(8), dimension(nzout)   , intent(in)        :: zout
    real(8), dimension(np,nzout), intent(out)       :: varout
    ! Locals 
    integer                                         :: i,k,k0
    real(8)                                         :: rat

    real, parameter     :: almost_zero = 1.e-6


    varout=null_value

    do i=1,np

      if (abs(zin(i,kz)) <= almost_zero .and. abs(zin(i,nzin)) <= almost_zero) cycle ! dry nodes have all zins =0.



      ! interpolate z
      if (zout(i) < zin(i,kz)) then ! zout below bathymetry
        ! do nothing
      elseif (zout(i) >= zin(i,nzin)) then
        !if (debug .and. i==debug_node) write(*,*) 'getting surface!'
        write(*,*) 'getting surface!'
        varout(i,i) = varin(i,nzin)
      else

        ! interp1d
        k0=-1! layer below
        do k=kz,nzin-1
          if (zin(i,k) <= zout(i) .and. zout(i) < zin(i,k+1)) then
            k0=k
            exit
          endif
        enddo !k


        rat = (zout(i)-zin(i,k0)) / (zin(i,k0+1)-zin(i,k0))
        varout(i,i) = varin(i,k0)*(1.-rat)+varin(i,k0+1)*rat

      endif


    enddo !i=1,np


  end subroutine
