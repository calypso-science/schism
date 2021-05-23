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
    integer                                         :: i,k,kout,k0
    real(8)                                         :: rat

    real, parameter     :: almost_zero = 1.e-6

    ! zout and zin should be negative

    !zin(kz:nzin)
    ! k=kz (bottom)  (=1 if only sigma, no zlevels below)
    ! k=nzin (top)

    !zout(1:nzout)
    ! k=1 (top)
    ! k=nzout (bottom)

    ! for a given time instant, if a requested level is above water varout == surface
    ! if below, varout is interpolated.

    !varout=null_value  !PROBLEM: if dry, varout will be null_value instead of 0.0
    ! problem because we need velocities to be zero when dry?

    varout=null_value

    do i=1,np

      if (abs(zin(i,kz)) <= almost_zero .and. abs(zin(i,nzin)) <= almost_zero) cycle ! dry nodes have all zins =0.

      do kout=1,nzout

        ! interpolate z
        if (zout(kout) < zin(i,kz)) then ! zout below bathymetry
          varout(i,kout) =-99999! do nothing
        elseif (zout(kout) >= zin(i,nzin)) then
          !if (debug .and. i==debug_node) write(*,*) 'getting surface!'
          write(*,*) 'getting surface!'
          varout(i,kout) = varin(i,nzin)
        else

          ! interp1d
          k0=-1! layer below
          do k=kz,nzin-1
            if (zin(i,k) <= zout(kout) .and. zout(kout) < zin(i,k+1)) then
              k0=k
              exit
            endif
          enddo !k

          if (k0<0) then
            ! should never enter here because below and above zin cases were comtemplated before
            write(*,*) "couldn't find any layers for zlevel ", zout(kout), 'at node ', i
            write(*,*) "zin = ", zin(i,:)
            stop 'interpolatez - selfe2netcdf - ERR02'
          endif

          rat = (zout(kout)-zin(i,k0)) / (zin(i,k0+1)-zin(i,k0))
          varout(i,kout) = varin(i,k0)*(1.-rat)+varin(i,k0+1)*rat

        endif

        !if (debug .and. i==debug_node) then
        ! write(*,*) 'kout, zout(kout), k0, zin(i,k0), zin(i,k0+1), varin(i,k0), varin(i,k0+1), varout(i,kout) = ', &
        !           kout, zout(kout), k0, zin(i,k0), zin(i,k0+1), varin(i,k0), varin(i,k0+1), varout(i,kout)
        !endif
  

      enddo !kout

    enddo !i=1,np


  end subroutine
