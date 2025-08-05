        program inputreader
        integer nfrs,n,step
        CHARACTER(len=32) :: arg
        open(unit=5,file="protname.txt",status='old')
        read(5,*)
        read(5,*)n
        read(5,*)nfrs
        close(5)
        CALL getarg(1,arg)
        read (arg, '(I10)') step
        WRITE(*,*)"STEP: ",step
        nfrs = INT(nfrs/step)
        write(*,*)n,nfrs
!       nfrs=500 !testing!
        call umatrix(n,nfrs,step)

	End program inputreader

        subroutine umatrix(n,nfrs,step)
        integer i,j,k,l,m,imin,jmin,a,nca,n,nfrs,step,counter
        real, dimension(nfrs) :: rx,ry,rz,lavmsq
        real, dimension(0:nfrs) :: dipcorr
        real, dimension(0:nfrs) :: avdipcorr
        DOUBLE PRECISION :: init
        character(32)protname
        character(16)aa,ii
        open(unit=5,file="protname.txt",status='old')
        read(5,'(A)')protname
        close(5)
        protname=adjustl(protname)
        rx=0.0
        ry=0.0
        rz=0.0
        lavmsq=0.0

        l = dummy
        i=l
        WRITE(ii,*)i
        ii=adjustl(ii)

        !read from trajectory
        open(unit=13,file='deltar_traj_'//TRIM(ii),status='old')
        read(13,*)rx 
        READ(13,*)ry 
        READ(13,*)rz
        close(13)

        !Introduce the stepping function
        rx = rx(::step)
        ry = ry(::step)
        rz = rz(::step)

	dipcorr=0.0
	avdipcorr=0.0
	!calculate tcf of bond autocorrelation
	do j=0,nfrs/4
	    if(mod(j,1000).eq.0)write(*,*)"tcf frame:",j
	    do k=1,nfrs-j
!and now the tcf
	      dipcorr(j)=dipcorr(j)+((rx(k)*rx(j+k)+ry(k)*ry(j+k)+rz(k)*rz(j+k)))
	    end do
	    dipcorr(j)=dipcorr(j)/real(nfrs-j)
	    !dipcorr(l,j)=dipcorr(l,j)/dipcorr(l,0)
	end do !end j time loop

	i=l
	  WRITE(ii,*)i
	  ii=adjustl(ii)
	  init = dipcorr(0)
	  OPEN(unit=24,file='tcf_deltar_'//trim(ii),status='unknown')
	  DO k=0,nfrs/4
	    WRITE(24,*)step*2*k,dipcorr(k)/init !dipcorr(i,0)
	  END DO
	  CLOSE(24)

	end subroutine
