! -*- f90 -*-

!f2py -m interpvert -h interpvert.pyf interpvert.f90
!Generate library
!f2py  -c interpvert.pyf  interpvert.f90

      subroutine interph(var,px,py,x0,y0,idx,idy,vout,np,nx0,ny0)
      integer np,ip,ii,jj,nx0,ny0
      real var(ny0,nx0),vout(np)
      real px(np),py(np),idx,idy
      real sdx,sdy,sdx1,sdy1,ppx,ppy,fac,tmpint,tmp

      do ip=1,np
        ppx=idx*(px(ip)-x0)+1.
        ppy=idy*(py(ip)-y0)+1.
!        print*,px(ip),py(ip),ppx,ppy
        tmpint=0.
        fac=0.
        sdx=mod(ppx,1.0)
        sdy=mod(ppy,1.0)
        sdx1=(1.-sdx)
        sdy1=(1.-sdy)
        ii=min(max(1,int(ppx)),nx0)
        jj=min(max(1,int(ppy)),ny0)
        if (abs(var(jj,ii)).lt.9999) then
          tmp=sdx1*sdy1
          tmpint=tmpint+tmp*var(jj,ii)
          fac=fac+tmp
        endif
        if (jj.lt.ny0.and.abs(var(jj+1,ii)).lt.9999) then
          tmp=sdx1*sdy
          tmpint=tmpint+tmp*var(jj+1,ii)
          fac=fac+tmp
        endif
        if (ii.lt.nx0.and.abs(var(jj,ii+1)).lt.9999) then
          tmp=sdx*sdy1
          tmpint=tmpint+tmp*var(jj,ii+1)
          fac=fac+tmp
        endif
        if (ii.lt.nx0.and.jj.lt.ny0.and.abs(var(jj+1,ii+1)).lt.9999) then
          tmp=sdx*sdy
          tmpint=tmpint+tmp*var(jj+1,ii+1)
          fac=fac+tmp
        endif
        if (fac.gt.0) then
          tmpint=tmpint/fac
        else
          tmpint=0.
        endif
        vout(ip)=tmpint
      enddo
      return
      end subroutine