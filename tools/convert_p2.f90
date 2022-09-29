program convert_p2
 implicit none
 integer :: i,j,m,nnodes,nbnodes,nelem,ndom,nb1,nedges,iter,pvalue
 real*8, allocatable :: x(:,:),xe(:,:)
 integer, allocatable :: tri(:,:),edge(:,:),trie(:,:)
 character*40 :: gridfile
 character*42 :: pgridfile
 !
 call getarg(1,gridfile)
 
 open(unit=1,file=gridfile,status='old')
 read(1,*) nnodes,nelem,pvalue
 allocate(x(2,nnodes))
 read(1,*) x
 allocate(tri(3,nelem))
 read(1,*) tri
 read(1,*) nbnodes
 close(1)
 !
 allocate(edge(7,3*nelem))
 edge=0
 !
 call findEdges(tri,3,nelem,edge,nedges)
 write(6,*) 'nedges=',nedges
 allocate(trie(3,nelem))
 allocate(xe(2,nedges))
 do i=1,nedges
    xe(:,i)=(x(:,edge(1,i))+x(:,edge(2,i)))*0.5
    trie(edge(5,i),edge(3,i))=i+nnodes
    if (edge(4,i) > 0) then
       trie(edge(6,i),edge(4,i))=i+nnodes
    endif
 enddo
 !
 pgridfile='p2_'//trim(adjustl(gridfile))
 open(unit=10,file=pgridfile,form='formatted')
 write(10,*) nnodes+nedges,nelem,2
 do i=1,nnodes
    write(10,*) x(:,i)
 enddo
 do i=1,nedges
    write(10,*) xe(:,i)
 enddo
 do i=1,nelem
    write(10,"(6(1x,I10))") tri(:,i),trie(:,i)
 enddo
 write(10,*) nbnodes
 do i=1,nbnodes
    write(10,"(2(1x,I10))") i,1
 enddo
 close(10)
 !
 deallocate(x,tri,trie,xe)
 !

601 format('TITLE ="',a40,'"')
602 format('VARIABLES ="X", "Y","')
603 format('ZONE T="VOL_MIXED",N=',i8,' E=',i8,' ET=TRIANGLE', &
         ' F=FEPOINT')
604 format(3(I7,1X))
 !
end program convert_p2
!
! find all edges in a given
! unstructured grid
!
subroutine findEdges(cellCon,nv,ncells,etmp,nedges)
!
implicit none
!
! subroutine arguments
!
integer, intent(in) :: nv
integer,intent(inout) :: ncells
integer, intent(inout) :: cellCon(nv,ncells)
integer, intent(inout) :: etmp(5,nv*ncells)
integer, intent(inout) :: nedges
!
! local variables
!
integer :: ii,i,j,jp1,maxedges,m
integer, allocatable :: iptr(:),iflag(:)
integer :: eloc(2)
integer :: nvert,nnodes
!
! begin
!
maxedges=ncells*4
!
nnodes=0
do i=1,ncells
 do m=1,nv
  nnodes=max(nnodes,cellCon(m,i))
 enddo
enddo
!
allocate(iptr(nnodes))
allocate(iflag(ncells))
!
etmp=0
!
iptr=0
!
do i=1,ncells
   nvert=nv
   if (nvert > 3) then
      if (cellCon(3,i)==cellCon(4,i)) nvert=3
   endif
   do j=1,nvert
      jp1=mod(j,nvert)+1
      eloc(1)=cellCon(j,i)
      eloc(2)=cellCon(jp1,i)
      call insert_edge(j,iptr,etmp,eloc,i,nedges,nnodes,maxedges)
   enddo
enddo
!
deallocate(iptr)
deallocate(iflag)
!
return
end subroutine findEdges

subroutine insert_edge(iedge,iptr,edge,eloc,cellIndex,nedges,nnodes,maxedges)
implicit none
!
! subroutine arguments
!
integer, intent(in) :: iedge,nnodes,maxedges,eloc(2),cellIndex
integer, intent(inout) :: nedges
integer, intent(inout) :: iptr(nnodes)
integer, intent(inout) :: edge(7,maxedges)
!
! local variables
!
integer :: e1(2),e2(2),ip,te
!
! begin
!
e1=eloc
if (e1(1).gt.e1(2)) then
   te=e1(1)
   e1(1)=e1(2)
   e1(2)=te
endif
!
ip=iptr(e1(1))
!
checkloop: do while(ip > 0)
   e2=edge(1:2,ip)
   if (e2(1).gt.e2(2)) then
      te=e2(1)
      e2(1)=e2(2)
      e2(2)=te
   endif
   if (sum(abs(e1-e2))==0) then
      edge(4,ip)=cellIndex
      edge(6,ip)=iedge
      return
   endif
   ip=edge(7,ip)
enddo checkloop
!
nedges=nedges+1
edge(1:2,nedges)=eloc
edge(3,nedges)=cellIndex
edge(5,nedges)=iedge
edge(7,nedges)=iptr(e1(1))
iptr(e1(1))=nedges
!
return
end subroutine insert_edge
