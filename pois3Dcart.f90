program mpiGather_noSub
use mpi
implicit none
integer,parameter::						nX=51,&
															nY=51,&
															nZ=51,&
															nDim=3
															
double precision,parameter::	xStart=0.0,&
															yStart=0.0,&
															zStart=0.0,&
															xEnd=1.0,&
															yEnd=1.0,&
															zEnd=1.0,&
															w=1.8d0

double precision,dimension(nX,nY,nZ)				::u,f,uOld
character(10)				:: filename
integer,dimension(nDim) ::bs,es,bn,en,siz,nP,myDim,myCoord,src,dest

logical :: myperiod,bottom,top,east,west,north,south
integer,allocatable,dimension(:) ::disp,sizz
integer,allocatable,dimension(:,:) ::tArr
integer ::sizDP, liney,planeyz,planexz,planexy,block,req
integer												::i,j,k,cnt,nproc, ierr, id, STATUS(mpi_status_size), flag, comm3d
double precision							::t1,t2,t3,tol,uAvg,uAvgOld,uAvgAll,uAvgOldAll,sumor, it1, it2, x(nX), y(nY), z(nZ), dx ,dy, dz
!=======================================================================
!==============================MPI======================================
CALL mpi_init(ierr)
!tStart=mpi_wtime()
call mpi_cart_create(mpi_comm_world,nDim,(/4,3,2/),(/.false.,.false.,.false./),.true.,comm3d,ierr)
call mpi_cart_get(comm3d,nDim,myDim,myperiod,mycoord,ierr)
CALL mpi_comm_rank(comm3d,id,ierr)
CALL mpi_comm_size(mpi_comm_world,nproc,ierr)
!call mpi_cart_coords(comm3d,5,2,mycoords,ierr)
do i=1,nDim
call mpi_cart_shift(comm3d,i-1,1,src(i),dest(i),ierr)
enddo
!call mpi_cart_shift(comm3d,1,1,srcy,desty,ierr)
!call mpi_cart_shift(comm3d,2,1,srcz,destz,ierr)
!write(*,22)id,mycoord(1),mycoord(2),mycoord(3), mydim(1), mydim(2), mydim(3)!, myperiod
!22 format(7(i3,1x))
!print*,id,mycoords(1),mycoords(2)
!print*,id,srcz,destz

ALLOCATE(disp(0:nproc-1),sizz(0:nproc-1))
ALLOCATE(tArr(0:maxval(myDim)-1,3))

east	=(mycoord(1)==myDim(1)-1)
west	=(mycoord(1)==0)
north	=(mycoord(2)==myDim(2)-1)
south	=(mycoord(2)==0)
top		=(mycoord(3)==myDim(3)-1)
bottom=(mycoord(3)==0)

nP=(/nX,nY,nZ/)

do k=1,nDim

	it1=nP(k)/myDim(k)
	it2=myDim(k)-mod(nP(k),myDim(k))

	tArr(0,1)=1
	do i=0,myDim(k)-1
		if(i==it2) it1=it1+1
		tArr(i,2)=tArr(i,1)+it1-1
		tArr(i,3)=tArr(i,2)-tArr(i,1)+1
		if(i==myDim(k)-1) exit
		tArr(i+1,1)=tArr(i,2)+1	
	enddo


	do i=0,myDim(k)-1
		if(i==mycoord(k)) then
			bs(k)=tArr(i,1)
			es(k)=tArr(i,2)
			siz(k)=tArr(i,3)
		endif
	enddo

	bn(k)=bs(k)
	en(k)=es(k)
	
	if((west .and. k==1) .or. (south 	.and. k==2) .or. (bottom .and. k==3)) bn(k)=bs(k)+1
	if((east .and. k==1) .or. (north 	.and. k==2) .or. (top 	 .and. k==3)) en(k)=es(k)-1
	
enddo

CALL mpi_type_vector (siz(2), 1, nX, mpi_double_precision,liney ,ierr)
CALL mpi_type_commit(liney,ierr)
CALL mpi_type_extent (mpi_double_precision,sizDP,ierr)
CALL mpi_type_hvector(siz(3), 1, nX*nY*sizDP,liney,planeyz,ierr)
CALL mpi_type_vector (siz(3), siz(1), nX*nY, mpi_double_precision,planexz ,ierr)
CALL mpi_type_vector (siz(2), siz(1), nX, mpi_double_precision,planexy ,ierr)

CALL mpi_type_hvector (siz(3),1,nX*nY*sizDP,planexy,block,ierr)

CALL mpi_type_commit(planeyz,ierr)
CALL mpi_type_commit(planexz,ierr)
CALL mpi_type_commit(planexy,ierr)
CALL mpi_type_commit(block,ierr)

!do i=1,nDim
!write(*,'(6(1x,i3))') id ,bs(i), es(i), bn(i),en(i),siz(i)
!CALL mpi_barrier(mpi_comm_world,ierr)
!if(id==0) write(*,*)
!CALL mpi_barrier(mpi_comm_world,ierr)
!enddo

!stop
!if(id==6) write(*,'(7(1x,i3))') id ,bn(1), en(1), bn(2),en(2),siz(1),siz(2)
CALL mpi_barrier(mpi_comm_world,ierr)
!stop
!-----------------------------MPI---------------------------------------

x(1)=xStart
dx = (xEnd-xStart)/(nX-1)
do i=1,nX-1
	x(i+1)=x(i)+ dx
enddo

y(1)=yStart
dy = (yEnd-yStart)/(nY-1)
do i=1,nY-1
	y(i+1)=y(i)+ dy
enddo

z(1)=zStart
dz = (zEnd-zStart)/(nZ-1)
do i=1,nZ-1
	z(i+1)=z(i)+ dz
enddo
!-----------------------------------------------------------------------
f(bs(1):es(1),bs(2):es(2),bs(3):es(3))=-5.0d0

u(bs(1):es(1),bs(2):es(2),bs(3):es(3))=0.0d0

t1=1.0/(x(2)-x(1))**2.0
t2=1.0/(y(2)-y(1))**2.0
t3=1.0/(z(2)-z(1))**2.0

tol=1e-3
flag=0
cnt=0
do while(.true.) !cnt .lt. 250
cnt=cnt+1
uOld(bs(1):es(1),bs(2):es(2),bs(3):es(3))=u(bs(1):es(1),bs(2):es(2),bs(3):es(3))

do k=bn(3),en(3)
do j=bn(2),en(2)
do i=bn(1),en(1)
u(i,j,k)=(1-w)*u(i,j,k) + w*( (u(i-1,j,k)+u(i+1,j,k))*t1 + (u(i,j-1,k)+u(i,j+1,k))*t2 + (u(i,j,k-1)+u(i,j,k+1))*t3-f(i,j,k) ) / (2*(t1+t2+t3))
!u(i,j)=( (u(i-1,j)+u(i+1,j))*t1+(u(i,j-1)+u(i,j+1))*t2-f(i,j) ) / (2*(t1+t2))
enddo
enddo
enddo

!send and receive
	call mpi_sendrecv(u(en(1),bn(2),bn(3)),1,planeyz,dest(1),50,u(bn(1)-1,bn(2),bn(3)),1,planeyz,src(1) ,50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(u(bn(1),bn(2),bn(3)),1,planeyz,src(1) ,50,u(en(1)+1,bn(2),bn(3)),1,planeyz,dest(1),50,mpi_comm_world,STATUS,ierr)
	
	call mpi_sendrecv(u(bn(1),en(2),bn(3)),1,planexz,dest(2),50,u(bn(1),bn(2)-1,bn(3)),1,planexz,src(2) ,50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(u(bn(1),bn(2),bn(3)),1,planexz,src(2) ,50,u(bn(1),en(2)+1,bn(3)),1,planexz,dest(2),50,mpi_comm_world,STATUS,ierr)
	
	call mpi_sendrecv(u(bn(1),bn(2),en(3)),1,planexy,dest(3),50,u(bn(1),bn(2),bn(3)-1),1,planexy,src(3) ,50,mpi_comm_world,STATUS,ierr)
	call mpi_sendrecv(u(bn(1),bn(2),bn(3)),1,planexy,src(3) ,50,u(bn(1),bn(2),en(3)+1),1,planexy,dest(3),50,mpi_comm_world,STATUS,ierr)

uAvg=0.0d0
uAvgOld=0.0d0
do k=bn(3),en(3)
do j=bn(2),en(2)
do i=bn(1),en(1)
	uAvg=uAvg+u(i,j,k)**2
	uAvgOld=uAvgOld+uOld(i,j,k)**2
enddo
enddo
enddo

CALL MPI_Reduce(uAvg   ,uAvgAll   ,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
CALL MPI_Reduce(uAvgOld,uAvgOldAll,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

if(id==0) then
	uAvgAll=sqrt(uAvgAll)/(nX*nY*nZ)
	uAvgOldAll=sqrt(uAvgOldAll)/(nX*nY*nZ)
	
	if(cnt==1) sumor=abs(uAvgAll-uAvgOldAll)

	if(abs(uAvgAll-uAvgOldAll)/sumor<tol) then !cnt>150
		print*,'Poisson exited at',cnt,'iteration : Residual',abs(uAvgAll-uAvgOldAll)
		flag=1
	endif

	if(abs(uAvgAll-uAvgOldAll)/sumor>1000) then !cnt>150
		print*,'Poisson diverged at',cnt,'iteration : Residual',abs(uAvgAll-uAvgOldAll)
		print*,'==============================================================='
		stop
	endif

	!print*,cnt,abs(uAvgAll-uAvgOldAll)/sumor
endif

call mpi_bcast(flag,1,mpi_integer,0,mpi_comm_world,ierr)
if(flag==1) exit

enddo
!=======================================================================
!print*,id,flag,cnt,abs(uAvgAll-uAvgOldAll)/sumor
!CALL mpi_barrier(mpi_comm_world,ierr)

!call mpi_gather(u(1,1,bs(id)),siz(id),oneplane,u(1,1,1),siz(id),oneplane,0,mpi_comm_world,ierr)
!call mpi_Allgather(u(1,1,bs(id)),siz(id),oneplane,u(1,1,1),siz(id),oneplane,mpi_comm_world,ierr)

!disp(0)=0
!do i=0,nproc
!if(i==id) sizz(i)=siz(3)
!if(i==nproc-1) exit
!disp(i+1)=disp(i) + sizz(i)
!enddo

!if(id==0) print*,'a'

!if(id/=0) then
!	if(mod(id,2)==0) then
!		CALL mpi_send(u(bs(1),bs(2),bs(3)),1,block,0,50,mpi_comm_world,ierr)
!		CALL mpi_recv(u(bs(1),bs(2),bs(3)),1,block,id,50,mpi_comm_world,STATUS,ierr)
!	else
!		CALL mpi_recv(u(bs(1),bs(2),bs(3)),1,block,id,50,mpi_comm_world,STATUS,ierr)
!		CALL mpi_send(u(bs(1),bs(2),bs(3)),1,block,0,50,mpi_comm_world,ierr)
!	endif
!endif

!print*,id			
!call mpi_Allgatherv(u(bs(1),bs(2),bs(3)),siz(3),planexy,u(1,1,1),sizz,disp,planexy,mpi_comm_world,ierr)
!=======================================================================
WRITE(filename,'(a,i2.2,a)')'3D',id+1,'.dat'
OPEN (UNIT=15,FILE=filename)

if(id==0) then
	write(15,*) 'Title = "Poisson Solution"'
	write(15,*) 'Variables = "x","y","z","u"'
	write(15,*) 'Zone k=',nZ,',j=', nY,',i=', nX,', DATAPACKING="POINT"'
endif

	do k=bn(3),en(3)
		write(15,*)
		do j=bn(2),en(2)
			write(15,*)
			do i=bn(1),en(1)
				write(15,101) x(i),y(j),z(k),u(i,j,k)
			end do
		end do
	enddo
	
	close(15)
	
	101 format(4(F20.12,1X))
	
CALL mpi_type_free(planeyz ,ierr)
CALL mpi_type_free(planexz,ierr)
CALL mpi_type_free(planexy,ierr)
CALL mpi_type_free(liney,ierr)
CALL mpi_type_free(block,ierr)
CALL mpi_finalize(ierr)

end program mpiGather_noSub
