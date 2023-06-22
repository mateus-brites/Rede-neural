c Main code to test subroutine sysredneu
c June 6, 2023

      implicit none
      integer nnx, nnobs, nrowsa, ndepth, nnAb,nnobs2
      parameter(nnx=100, nnobs=2000, nrowsa=100, ndepth=10, nnAb=10000)
      parameter(nnobs2 = nnobs*2)
      integer nobs
      integer nx, depth, ncomps(ndepth), npsi
      integer nAb
      integer ncoan
      double precision Ab(nnAb), y(nnobs), f(nnobs),vaux(nnobs2)
      double precision jacob(nnobs, nnAb)
      double precision jacosave(nnobs,nnAb), fsave(nnobs)      
      double precision jaca(nrowsa,nnAb), jacx(nnAb,nnAb)
      double precision aat(nnobs,nnobs), aatsave(nnobs,nnobs)
      double precision jacaux(nrowsa,nnAb)
      double precision a(nrowsa, nrowsa, ndepth), b(nrowsa, ndepth)
      double precision x(nnx, nnobs)
      double precision seed
      integer iperm(nnobs)
      integer i, j, k, iobs, idep, kprob
      integer kon, max, nvar, nref
      integer idata
      integer iadd
      real day1, day2
      double precision eps 

      double precision z, under
      integer sorte, iqn, maxrow, ipri, inewran, option

c  Data that will be read from fort.106
      double precision datos(4,6060)

 
      seed = 28383723771717.d0
      write(*, *)' Number of function psi:'
      read(*, *) npsi

      write(*, *)' If want detailed printing type 1, else type 0'
      read(*, *) ipri
      


      write(*, *)' Maximal number of observations handled at each
     * iteration:'
      read(*, *) maxrow
      write(*, *)' Maximal number of observations handled at each
     * iteration:',maxrow    
 

      write(*, *)' Under-relaxation parameter (default:1):'
      read(*, *) under
      write(*, *)' Under-relaxation parameter:', under

      write(*, *)' Number refinements (recall that 0 is "pure newton"):'
      read(*, *) nref
      write(*, *)' Refinements following each Newton iteration:', nref

      if(nref.gt.1) then
      write(*, *)' If want quasi-Newton refinements type 1, else type 0'
      read(*, *) iqn
      endif


      write(*, *)' Maximal number of iterations:'
      read(*, *) max

      write(*, *)' Stopping epsilon:'
      read(*, *) eps
  
 

      write(*, *)' Data generation:'
      write(*, *)' If you want Saint Venant given in fort.106, type 1'
      write(*, *)' Otherwise, type 2'
      read(*, *) kprob
      if(kprob.eq.1) then

c  Saint Venant

      nvar = 3
     
      do iobs = 1, 6060
      read(106, *) (datos(i, iobs),i=1,4)
      end do




      write(*, *)' Observations given every 12 hours'    
      write(*, *)' Number of observations nobs:'

      read(*, *) nobs       
      if(nx.gt.nnx.or.nobs.gt.nnobs) then
      write(*, *)' Abort. Maximal nx and nobs are:', nnx, nnobs
      stop
      endif
      if(2*nobs.gt.nnobs) then
      write(*, *)' Aborted because 2*nobs must be < nnobs'
      stop
      endif 


      write(*, *)' Training data will go from day: '
      read(*, *) day1
      write(*, *)' to day: '
      read(*, *) day2
      write(*, *)' Training data are between days', day1,' and', day2
      do iobs = 1, nobs
3     call randin(seed, 1, 6060, j)
      if((datos(3,j).lt.day1).or.(datos(3,j).gt.day2))
     *   go to 3
      x(1, iobs) = datos(1, j)
      x(2, iobs) = datos(2, j)
      x(3, iobs) = datos(3, j)
      y(iobs) = datos(4, j)     
      end do
      write(*, *)' Test data will go from day: '
      read(*, *) day1
      write(*, *)' to day: '
      read(*, *) day2

      do iobs = nobs+1, 2*nobs
5     call randin(seed, 1, 6060,j)
      if((datos(3,j).lt.day1).or.(datos(3,j).gt.day2))
     *   go to 5        
      x(1, iobs) = datos(1, j)
      x(2, iobs) = datos(2, j)
      x(3, iobs) = datos(3, j)
      y(iobs) = datos(4, j)
      end do 

      write(*, *)' Test data are between days', day1,' and', day2    

      
c      do iobs = 1, 2*nobs
c      call randin(seed, 1, 6060, j) 
c      x(1, iobs) = datos(1, j)
c      x(2, iobs) = datos(2, j)
c      x(3, iobs) = datos(3, j)
c      y(iobs) = datos(4, j)
c      end do

      write(*, *)' Data for Saint Venant (Qinlet, x, days, z):'
      do iobs = 1, 2*nobs
      write(*, *) iobs, (x(i, iobs),i=1,3), y(iobs)
      end do



      





      write(*, *)' Observation data have been read'



      write(*, *)' If you want to add quadratic dependence type 1'
      write(*, *)' Otherwise type 2'
      read(*, *) i
      if(i.ne.1) then
      nx = nvar
      else
      write(*, *)' Variables x include quadratic terms.'
      do iobs=1, 2*nobs
      x(nvar+1, iobs) = x(1, iobs)**2
      x(nvar+2, iobs) = x(2, iobs)**2
      x(nvar+3, iobs) = x(3, iobs)**3

      x(nvar+4, iobs) = x(1, iobs)**2
      x(nvar+5, iobs) = x(2, iobs)**2
      x(nvar+6, iobs) = x(3, iobs)**2 
      x(nvar+7, iobs) = x(1, iobs)**2 * x(2, iobs)
      x(nvar+8, iobs) = x(1, iobs)**2 * x(3, iobs)
      x(nvar+9, iobs) = x(2, iobs)**2 * x(3, iobs)
      x(nvar+10, iobs) = x(1, iobs) * x(2, iobs)**2
      x(nvar+11, iobs) = x(1, iobs) * x(3, iobs)**2
      x(nvar+12, iobs) = x(2, iobs) * x(3, iobs)**2      
      x(nvar+13, iobs) = dsin(x(1, iobs))
      x(nvar+14, iobs) = dsin(x(2, iobs))
      x(nvar+15, iobs) = dcos(x(3, iobs))
      x(nvar+16, iobs) = dcos(x(1, iobs))
      x(nvar+17, iobs) = dcos(x(2, iobs))
      x(nvar+18, iobs) = dcos(x(3, iobs))
 



      end do

      write(*, *)' How many additional input xs you want?'
      write(*, *)' The first one is Q-inlet squared'
      read(*, *) iadd

      nx = nvar+iadd
      write(*, *)' nx including quadratic dependences:', nx
      endif 
             


      write(*, *)' Reading Saint Venant data is complete'
      endif
     
      if(kprob.ne.1) then
      write(*, *)' Dimension of x :'
      read(*, *) nx
      write(*, *)' Number of observations nobs:'
      read(*, *) nobs
      if(nx.gt.nnx.or.nobs.gt.nnobs) then
      write(*, *)' Abort. Maximal nx and nobs are:', nnx, nnobs
      stop
      endif

      write(*, *)' nx = ', nx,' nobs = ', nobs

      write(*, *)' Data generation'
      write(*, *)' Type 1 for all-random generation'
      write(*, *)' Type 2 for function ||x||^2'
      write(*, *)' Type 3 for function sum sin(x_i)'
      read(*, *) idata
      if(idata.eq.1) then
c  All random generation
c  We generate 2 nobs observations in order to use
c  the last nobs as Test Set

      write(*, *)' All-random data generation'      
      do iobs = 1, nobs+nobs
      do i = 1, nx
      call rondo(seed, x(i,iobs))
      end do
      call rondo(seed, y(iobs))
      end do
      endif

      if(idata.eq.2.or.idata.eq.3) then
      do iobs = 1, nobs+nobs
      do i = 1, nx
      call rondo(seed, x(i,iobs))
      end do
      y(iobs) = 0.d0
      do i = 1, nx
      if(idata.eq.2) then
      y(iobs) = y(iobs) + x(i, iobs)**2
      endif
      if(idata.eq.3) then
      y(iobs) = y(iobs) + dsin(x(i, iobs))
      endif
      end do
      end do
      endif
      endif
c  This endif corresponds to "if(kprob.ne.1)"



      write(*, *) ' Depth = '
      read(*, *) depth
      if(depth.gt.ndepth) then
      write(*, *)' Abort. Maximal depth is:', ndepth
      stop
      endif

      write(*, *)' depth = ', depth

      do i = 1, depth
      write(*, *)' Components of level', i,'  (Rows of A number', i,'):'
      if(i.eq.depth) then
      ncomps(i) = 1
      else
      read(*, *) ncomps(i)
      endif
      if(ncomps(i).gt.nrowsa) then
      write(*, *)' Abort, Maximal number of components at each level:',
     *  nrowsa
      stop
      endif
      end do
      ncomps(depth) = 1
      write(*, *)
      do i = 1,depth
      write(*, *)' Components of level', i,'  (Rows of A number', i,'):'
     *  , ncomps(i)   
      end do
      

      nAb = (nx+1)*ncomps(1)
      do i = 2, depth
      nAb = nAb + (ncomps(i-1)+1)*ncomps(i)
      end do



      write(*, *)' Total number of coefficients (A, b) :', nAb
      if(nAb.gt.nnAb) then
      write(*, *)' Abort. Total number of coeffs (A, b) exceeds', nnAb
      stop
      endif

      write(*, *)' To continue type 1'
      read(*, *) i
 


      write(*, *)
      write(*, *) ' Data:'
      write(*, *)' iobs   x_1,...x_{nx}   y'
      do iobs = 1, nobs
      write(*, *) iobs, (x(i, iobs),i=1,nx), y(iobs)
      end do

      write(*, *)

6     write(*, *)' To continue type 1'
      read(*, *) i

      write(*, *)' Generation of Initial NN for sysredneu'
      write(*, *)' If you want a random generation type 1'
      write(*, *)' If you want to use the NN that comes '
      write(*, *)' from previous run, type 2'

      read(*, *) inewran
      if(inewran.eq.2) then
      write(*, *) ' Initial NN comes from previous run'
      else
      write(*, *)'  Initial NN is random'
      endif 
 

c Generation of coefficients (A, b)
      if(inewran.eq.1) then
      ncoan = nx

      do k = 1, depth
      do i = 1, ncomps(k)
      do j = 1, ncoan 
      call rondo(seed, a(i, j, k))
      end do
      call rondo(seed, b(i, k))
      b(i, k) = b(i, k)+1.d0
      end do
      ncoan = ncomps(k)
      end do
      endif

      if(inewran.eq.2) then
      write(*, *) ' The new run adds a coordinate to x'
      nx = nx+1
      do i = 1, ncomps(1)
      a(i, nx, 1) = 0.d0
      end do
      write(*, *)' Decide which should be the new coordinate of x'
      write(*, *)' Option 1:', ' Previous x(1) squared'
      write(*, *)' Option 2:', ' Previous x(2) squared'
      write(*, *)' Option 3:', ' Product x(1)*x(2)'   
      read(*, *) option
      if(option.eq.1) then
      do iobs = 1, 2*nobs
      x(nx, iobs) = x(1, iobs)**2
      end do
      endif
      if(option.eq.2) then
      do iobs = 1, 2*nobs
      x(nx, iobs) = x(2, iobs)**2
      end do
      endif
      if(option.eq.3) then
      do iobs = 1, 2*nobs
      x(nx, iobs) = x(1, iobs) * x(2, iobs)
      end do
      endif       
      endif
     
      if(ipri.eq.1) then
      ncoan = nx   
      do k = 1, depth
      
      write(*, *)' Matrix A number', k
      do i = 1, ncomps(k)
      write(*, *)(a(i, j, k),j=1,ncoan)
      end do
      ncoan = ncomps(k)
      end do
      
      do k = 1, depth
      write(*, *)' Vector b number', k
      write(*, *)(b(i, k),i=1,ncomps(k))
      end do
      endif

      k = 1
      ncoan = nx
      do idep = 1, depth
      do i = 1, ncomps(idep)
      do j = 1, ncoan
      Ab(k) = a(i, j, idep) 
      k = k + 1
      end do
      Ab(k) = b(i, idep) 
      k = k + 1
      end do
      ncoan = ncomps(idep)
      end do
      if(ipri.eq.1) then
      write(*, *)' Vector Ab:'
      write(*, *)(Ab(i),i=1,nAb)
      endif

      maxrow = min0(maxrow, nobs)

2     call sysredneu(x,nx,y,nobs,depth,ncomps,nAb,Ab,npsi,f,jacob,
     *jacosave,fsave,nnobs,nnx,nrowsa,jaca,jacx,jacaux,a,b,kon,max,seed,
     *under,eps,maxrow,nref,iqn,ipri,aat,aatsave,vaux,iperm)

      write(*, *)' Number of iterations performed:', kon

      write(*, *)' Recall that:'

      write(*, *)' nx = ', nx,' nobs = ', nobs
      write(*, *)' depth = ', depth
      write(*, *)' Number of function psi:', npsi

      write(*, *)' Total number of coefficients (A, b) :', nAb         

      do i = 1,depth
      write(*, *)' Components of level', i,'  (Rows of A number', i,'):'
     *  , ncomps(i)   
      end do
     
c  Transform the vector Ab into the matrices a(i,j,idep) and b(i,idep) 
      k = 1
      ncoan = nx
      do idep = 1, depth
      do i = 1, ncomps(idep)
      do j = 1, ncoan
      a(i, j, idep) = Ab(k) 
      k = k + 1
      end do
      b(i, idep) = Ab(k) 
      k = k + 1
      end do
      ncoan = ncomps(idep)
      end do

      write(*, *)' If you want to expand the vector x using the'
      write(*, *)' the actual NN as initial point, type 1'
      read(*, *) i
      if(i.eq.1) then
      write(*, *)' We use the NN output of sysredneu to initialize'
      write(*, *)' new optimization process'
      go to 6
      endif
 




      if(max.eq.0) go to 4  

      write(*, *)' To continue with Conference Test, type 1'
      read(*, *) i
      if(i.ne.1) stop

      maxrow = nobs


      max = 0
      go to 2

4     write(*, *)' To continue with Test Set, type 1'
      read(*, *) i
      if(i.ne.1) stop

      maxrow = nobs
      do iobs = 1, nobs
      do i = 1, nx
      x(i, iobs) = x(i, nobs+iobs)
      end do
      y(iobs) = y(nobs+iobs)
      end do

      max = 0
      go to 2          


      stop
      end


 


      subroutine sysredneu(x,nx,y,nobs,depth,ncomps,nAb,Ab,npsi,f,jacob,
     *jacosave,fsave,nnobs,nnx,nrowsa,jaca,jacx,jacaux,a,b,kon,max,seed,
     *under,eps,maxrow,nref,iqn,ipri,aat,aatsave,aux,iperm)
      implicit none
      integer nnobs, nrowsa, nnx 
      integer nx, npsi, nref
      integer nobs
c  positions for aux must be at least 2*nobs      
      double precision x(nnx, nobs), y(nobs), aux(2*nobs)
      double precision jaca(nrowsa, *),jacx(nrowsa, *),jacaux(nrowsa,*)
      double precision aat(nnobs, nobs), aatsave(nnobs,nobs), add
      integer iperm(nobs)
      logical posdef
      
      double precision ff(1)

      integer ndats, nAb, one
      double precision Ab(nAb), f(nobs), jacob(nnobs,nAb)
      double precision jacosave(nnobs,nAb), fsave(nobs)
      integer k, idep, depth, ncoan
      integer ncomps(depth)
      integer i, j
      parameter(one=1)
      double precision a(nrowsa, nrowsa, depth), b(nrowsa, depth)
      integer iobs, kon
      integer modo
      double precision under
      integer max, mlin, iermgs
      integer sorte
      integer maxrow
      integer iqn, ipri


      double precision s(nAb)
      double precision snor, z, seed

      double precision error, err2, rmsd, eps
      double precision sts
 
      integer discar

c The observations are:
c     x(1, 1), ..., x(nx, 1),    y(1)
c     x(1, 2), ..., x(nx, 2),    y(2)
c ................................
c     x(1, nobs)...,x(nx, nobs), y(nobs)
c
c (This table explains the meaning of the inputs x, y, nx, and nobs.
c
c  Attention: Notice that each component of y is a single number. 
c  This means that we are restricted to Neural Networks that are scalar functions
c  In our notation ncomps(depth) is always equal to 1. 


c nobs is the number of observations
c nAb is the number of parameters (A, b) of the Neural Network
c Ab is the vector of parameters (A, b) in the form 
c ( row corresponding to A_1, b_1, ..., row corresponding to A_depth b_depth)
c where each "row corresponding to A_k, b_k 
c is is displayed in the form (row 1 of A_k, entry 1 of b_k, ..., row ncomp(k) of A_k, ..., entry ncomp(k) of b_k)
c 

      do iobs = 1, nobs
      iperm(iobs) = iobs
      end do

      if(maxrow.gt.nobs) maxrow=nobs


      if(ncomps(depth).ne.1) then
      write(*, *)' Error.In this implementation ncomps(depth) must be 1'
      stop
      endif

      kon = 0


c Evaluation of the system, obtaining F(Ab) and Jacobian(Ab)   
c Transformation of vector Ab into matrices A and vectors b

      write(*, *)' In this problem nx = ', nx,'  nobs=', nobs
      write(*, *)' Number of coefficients to be estimated:', nAb

      if(nAb.lt.nobs) then
      write(*, *)' nAb = ', nAb,' nobs =', nobs
      write(*, *)' Warning! This method was firstly developed ',
     * ' for nAb .ge. nobs!'
      endif
 
      write(*, *)' To continue type 1'
      read(*, *) i

     


c  Transform the vector Ab into the matrices a and the vectors b
1     k = 1
      ncoan = nx
      do idep = 1, depth
      do i = 1, ncomps(idep)
      do j = 1, ncoan
      a(i, j, idep) = Ab(k)
      k = k + 1
      end do
      b(i, idep) = Ab(k)
      k = k + 1
      end do
      ncoan = ncomps(idep)
      end do

      write(*, *)
      write(*, *)' Iteration:', kon

c      write(*, *)' Ab=', (Ab(i),i=1,nAb)

      write(*, *)' Ab(1)=', Ab(1),' ... Ab(',nAb,') =', Ab(nAb)
      write(*, *)



c Effective evaluation of the system , obtaining F(Ab) and Jacobian(Ab)

      if(mod(kon, nref+1).eq.0) then
      modo = 2
      else
      modo = 1
      endif

      if(maxrow.lt.nobs) then
      do iobs = 1, maxrow
      call randin(seed,iobs,nobs,sorte)
      i = iperm(iobs)
      iperm(iobs) = iperm(sorte)
      iperm(sorte) = i
      end do      
      endif
 


c  System and perhaps Jacobian evaluation

      do iobs = 1, maxrow 

      
      call redneu(a, b, x(1, iperm(iobs)), nx, one, npsi, depth, 
     *   ncomps,ff(1), jaca, jacx, jacaux, modo)

      

      if(ipri.eq.1) then
      write(*, *)' Observation ', iperm(iobs),' Predicted:', ff(1),
     *     ' Observed:', y(iperm(iobs))
      endif

      f(iobs) = ff(1) - y(iperm(iobs))

      if(modo.eq.2) then
      do j = 1, nAb
      jacob(iobs, j) = jaca(ncomps(depth), j)
      end do
      endif

      end do

c End of system and (perhaps)  Jacobian evaluation

      if(maxrow.ge.nobs.and.modo.eq.1.and.iqn.eq.1) then
c Proceed to update the Jacobian using Broyden quasi-Newton
c Previous Jacobian is in jacosave (nobs rows, nAb columns)
c Previous F is in fsave (nobs components)
c Last increment is in vector s (nAb components)
      sts = 0.d0
      do j = 1, nAb
      sts = sts + s(j)*s(j)
      end do
      do iobs = 1, maxrow
      aux(iobs) = f(iobs)-fsave(iobs) 
      do j = 1, nAb
      aux(iobs) = aux(iobs) -  jacosave(iobs, j)*s(j)
      end do
      aux(iobs) = aux(iobs)/sts
      end do
      do j = 1, nAb
      do iobs = 1, maxrow
      jacob(iobs, j) = jacosave(iobs, j) + aux(iobs)*s(j)
      end do
      end do
      endif


c  Recall that ncomps(depth) = 1 in this implementation

      do iobs = 1, maxrow 
      if(ipri.eq.1) then
      write(*, *)' f(', iperm(iobs),') = ', f(iobs),
     *   ' for data =', y(iperm(iobs))
      endif
      end do



c  Test stopping criterion
      if(maxrow.eq.nobs) then
      error = 0.d0
      err2 = 0.d0
      do iobs = 1, nobs
      error = dmax1(error, dabs(f(iobs)))
      err2 = err2 + f(iobs)**2
      end do

      rmsd = dsqrt(err2/dfloat(nobs))

      write(*, *)' Error (max diff between predicted & observed):',error
      write(*, *)' Sum of squares:', err2,' RMSD =', rmsd
      endif
      
      if(max.eq.0) then
      write(*, *)' Above resuts have been obtained for conference set'
      write(*, *)' or for test set, as maximum of iterations was 0'
      return
      endif  

      write(*, *)' Number of iterations up to now:', kon             
 
      if(maxrow.eq.nobs) then
      if(error.le.eps) then
      write(*, *)' Return by small error'
      return
      endif
      endif

      if(kon.ge.max) then
      write(*, *)' Number of iterations up to now:', kon     
      write(*, *)' Return by maximum of iterations', max,' reached'
      return
      endif 

c***********************************************************************
      if(nobs.le.nAb) then

c      write(*, *)' Parameters entry of mgs:'
c      write(*, *)' nobs=', nobs,' nAb=', nAb

      do iobs = 1, maxrow 
      do j = 1, nAb
      jacosave(iobs, j) = jacob(iobs, j)
      end do
      fsave(iobs) = f(iobs)
      end do

c      write(*, *)' Independent term that enters to mgs:',(f(i),i=1,maxrow)
c      write(*, *)' Matrix jacob that enters to mgs:'
      
c      do i = 1, nobs 
c      write(*, *)(jacob(i,j),j=1,nAb)
c      end do 



      
      z = 0.d0
      do j = 1, nAb
      do i = 1, maxrow
      z = dmax1(z, dabs(jacob(i,j)))
      end do
      end do
      write(*, *)' Maximal element of ', maxrow,' rows of Jacob:', z

      call mgs(jacob, f, s, nobs, nAb, iermgs, nnobs, discar,maxrow)
      write(*, *)' Discarded rows at mgs:', discar,' Rank of Jacobian:', 
     *   maxrow-discar
      write(*, *)' Error code ier of mgs:', iermgs
 

      if(iermgs.ne.0) then
      do iobs = 1, nobs
      do j = 1, nAb
      jacob(iobs, j) = jacosave(iobs, j)
      end do
      f(iobs) = fsave(iobs)
      end do

      do j = 1, maxrow
      do i = 1, maxrow
      aat(i,j)=0.d0
        do k = 1, nAb
        aat(i,j)=aat(i,j) + jacosave(i,k)*jacosave(j,k)
        end do
        aatsave(i,j) = aat(i,j)
      end do
      end do

      z = 0.d0
      do i = 1, maxrow
      z = dmax1(z, dabs(aat(i,i)))
      end do

      if(z.eq.0.d0) then
      write(*, *)' Diagonal of AA^t is null' 
      z = 1.d0
      else
      z = 1.d-3*z
      endif
3     do i = 1, maxrow
      aat(i,i)=aatsave(i,i) + z
      end do 

      call chole (maxrow, aat, posdef, add, nnobs)    
      if(.not.posdef) then
      write(*, *)' For z = ', z,' Matrix aat is not posdef yet'
      z = 2.d0*z 
      goto 3 
      endif
      call sicho (maxrow, aat, aux, f, aux(nobs+1), nnobs) 

      do j = 1, nAb
      s(j) = 0.d0
      do iobs=1,maxrow
      s(j) = s(j)+jacob(iobs,j)*aux(iobs)
      end do
      end do  
      endif
c  This endif corresponds to "if(iermgs.ne.0)"
 

      do iobs = 1, nobs
      do j = 1, nAb
      jacob(iobs, j) = jacosave(iobs, j)
      end do
      f(iobs) = fsave(iobs)
      end do          

      endif
c  This endif corresponds to "if(nobs.le.nAb)"
c*************************************************************************
      if(nobs.gt.nAb) then
      write(*, *)' The case nobs.gt.nAb is under revision'
      stop
      endif
c*************************************************************************


      snor = 0.d0
      do i = 1, nAb
      snor = dmax1(snor, dabs(s(i)))
      end do
      write(*, *)' Sup norm of the increment s:', snor

      if(snor.le.eps*dfloat(maxrow)/dfloat(nobs)) then
      write(*, *)' Finish because small increment'
      write(*, *)' Iterations performed:', kon
      return
      endif

c  Save current iterate Ab  

      do i = 1, nAb
      s(i) = - under * s(i)
      Ab(i) = Ab(i) +  s(i)
      end do

 
      kon = kon + 1
      go to 1
      end
      

   	subroutine mgs (a, b, x, m, n, ier, mlin, discar,maxrow)
	implicit none 
      integer m, n, ier, mlin
      double precision a(mlin, n), b(m),  x(n)
      integer ii, i, j
      double precision anor, z, tol, tolanor
      double precision pesca
      double precision seed
      integer discar, maxrow
      double precision snor

c  Modified  Gram-Schmidt subroutine for solving A x = b
c  with A being m x n and m <= n.
c  The matrix and the independent term are  destroyed.
c  Computes minimum norm solution when the system is 
c  compatible. 
c  If the system is not compatible
c  the equations whose rows are dependent with respect
c  to other ones are eliminated.

      discar = 0

	ier = 0
	if(m.eq.0.or.n.eq.0) then
	ier = 1
	return
	endif

c  Compute the maximal norm of the rows of A
	
	anor = 0.d0
	do j=1,n
	z = 0.d0
	do i=1,maxrow  
	z = z + a(i,j)*a(i,j)
	end do
	z = dsqrt(z)
	anor = dmax1(anor, z)
	end do

	if(anor.eq.0.d0) then
      write(*, *)' Matrix entering mgs is null'
      discar = maxrow
	ier = 2
      do j = 1, n
      x(j) = 0.d0
      end do
	return
	endif

	tol = 1.d-12

	tolanor = tol*anor

	do i = 1, maxrow
c  Normalize row i
	z = 0.d0
	do j = 1, n
	z = z + a(i, j)*a(i, j)
	end do
      z = dsqrt(z)

	if(z.lt.tolanor) then
c  Row i of A is a linear combination of rows 1,...,i-1
c  We replace this row by (0, ..., 0) and we set, 
c  accordingly, b(i) equal to 0.
c  Warning is represented by ier = 3
      z = 0.d0
      discar = discar+1
      ier = 3
      do j = 1, n 
      a(i, j) = 0.d0
      end do
      b(i) = 0.d0
      else 
	
	do j = 1, n
	a(i, j) = a(i, j)/z
	end do
	
c   Same multiplication in independent term
	b(i) = b(i)/z
      endif
	

	if(i.eq.maxrow) go to 1

      if(z.ne.0.d0) then
	do ii = i+1,  m
c   Orthogonalize row ii with respect to row i
	pesca = 0.d0
	do j =1, n
	pesca = pesca + a(ii, j) * a(i, j)
	end do
	do j = 1, n
	a(ii, j) = a(ii, j) - pesca * a(i, j)
	end do
c   Same in independent term
	b(ii) = b(ii) - pesca * b(i)
	end do
      endif

	end do


c   Compute solution = Q(trasp) b
1	do j = 1, n
	x(j)= 0.d0
	do i = 1, maxrow
	x(j) = x(j) + a(i, j) * b(i)
	end do
	end do


	return
	end



 
 



      subroutine redneu(a, b, x, n, m, npsi, depth, ncomps, f, jaca,
     *   jacx, jacaux, modo) 

      implicit none
      integer ndim  
      parameter(ndim = 100)

      double precision a(ndim,ndim,*), b(ndim,*),
     *    jaca(ndim, *), jacx(ndim,*), jacaux(ndim,*)
      double precision x(n), f(m)
      integer depth
      integer ncomps(depth)
      double precision y(ndim)
      integer mmax, nmax, i, nivel, j, k  
      integer n, m
      integer ncolult, ncpan
      double precision z(ndim,ndim), ynew(ndim)
      integer nada
      integer npsi
      integer modo
      double precision dpsi


      if(modo.eq.1) then

      do i = 1, n
      y(i) = x(i)
      end do
      ncpan = n

c  Compute New y
      do nivel = 1, depth
      do i = 1, ncomps(nivel)
      ynew(i) = b(i, nivel)
      do j = 1, ncpan
      ynew(i) = ynew(i) + a(i, j, nivel)*y(j)
      end do
      end do
      ncpan = ncomps(nivel)
      do i = 1, ncomps(nivel)
      call redneupsi(npsi, ynew(i), y(i), dpsi, nivel, depth) 
      end do           
      end do 




      do i = 1, ncomps(depth)
      f(i) = y(i)
      end do

      return

      endif
c  This endif corresponds to if(modo.eq.1)
c  From now on we have modo=2, so we compute both F and Jacobian
c 


      mmax = 0
      do i = 1, depth
      mmax = max0(mmax, ncomps(i))
      end do

      nmax = n*ncomps(1)
      do i = 2, depth
      nmax = nmax + ncomps(i-1)*ncomps(i) 
      end do

      do j = 1, nmax
      do i = 1, mmax
      jaca(i, j) = 0.d0
      jacaux(i, j) = 0.d0
      end do
      end do

      do j = 1, n 
      do i = 1, n
      jacx(i, j) = 0.d0
      end do
      jacx(j, j) = 1.d0
      end do




      do nivel = 1, depth 



c      write(*, *) ' Update Jaca at level ', nivel 

c      read(*, *) nada   



      if(nivel.eq.1) then

      do i = 1, n
      y(i) = x(i)
      end do

      ncpan = n
      ncolult = 0
      else
      ncpan = ncomps(nivel-1)
      endif

c  Use chain rule in order to update pre-existing Jacobian
c------------------------------------------------------------------

       
      if(nivel.gt.1) then


c  Columns from j=1 to j = ncolult must be updated

c      write(*, *)' Level ',nivel,' Update cols 1 to ',ncolult,' of jaca'

c      read (*, *) nada

c      write(*, *)' Columns 1 to ', ncolult, ' of jaca before update:'
c      do i = 1, ncpan
c      write(*, *) (jaca(i, j), j = 1, ncolult)
c      end do 

c      read(*, *) nada

c      write(*, *)' Columns 1 to ncomps(nivel-1) of A:'
c      do i = 1, ncomps(nivel)
c      write(*, *) (a(i, j, nivel),j=1,ncpan)
c      end do

c      read (*, *) nada


      do j = 1, ncolult
      do i = 1, ncomps(nivel)
      jacaux(i, j) = 0.d0
      do k = 1, ncpan
      jacaux(i, j) = jacaux(i, j) + a(i, k, nivel)*jaca(k, j)
      end do
      end do
      end do 

      do j = 1, ncolult
      do i = 1, ncomps(nivel)
      jaca(i, j) = jacaux(i, j)
      end do
      end do 
c      write(*, *)
c      write(*, *)' Columns 1 to ', ncolult, ' of jaca after update:'
c      do i = 1, ncomps(nivel)
c      write(*, *) (jaca(i, j), j = 1, ncolult)
c      end do 

c      read(*, *) nada
 


      endif 
c  End of "if (nivel .ne. 1)"
c-------------------------------------------------------------------

c-------------------------------------------------------------------
c  Add the columns of jaca that correspond to the new level

c      write(*, *)' Add the new columns of jaca from col', ncolult+1

c      read(*, *) nada

c      write(*, *) ' nivel =', nivel,' ncomps(nivel)=', ncomps(nivel)
c      write(*, *) ' ncpan = ', ncpan

      do k = 1, ncomps(nivel)
      do j = 1, (ncpan+1)*ncomps(nivel)
      jaca(k, ncolult + j)=0.d0
      end do
      end do


      do k = 1, ncomps(nivel)
      do j = 1, ncpan + 1
      if(j.lt.ncpan+1) then
      jaca(k , ncolult + (k-1)*(ncpan+1) + j ) = y(j)
      else
      jaca(k,  ncolult + (k-1)*(ncpan+1) + j) = 1.d0
      endif
      end do
      end do


c      write(*, *) ' First column typed here is number: ', ncolult+1
c      do k = 1, ncomps(nivel)
c      write(*, *) (jaca(k, j), j = ncolult+1, ncolult +
c     *     (ncomps(nivel)-1)*(ncpan+1)+ncpan+1) 
c      end do

c--------------------------------------------------------------------
c      read(*, *) nada 

      k = ncomps(nivel)
c      write(*, *)' Old ncolult:', ncolult
      ncolult =  ncolult + (k-1)*(ncpan+1)+ncpan+1
c      write(*, *)' New ncolult:', ncolult

c      read(*, *) nada

   
c--------------------------------------------------------------------
c  Update Jacobian with respect to x
      do j = 1, n
      do i = 1, ncomps(nivel)
      z(i, j) = 0.d0
      do k = 1, ncpan
      z(i, j) = z(i, j) + a(i, k, nivel)*jacx(k, j)
      end do
      end do
      end do

      do j = 1, n
      do i = 1, ncomps(nivel)
      jacx(i, j) = z(i, j)
      end do
      end do

c      write(*, *)' Updated Jacobian wrt x at level ', nivel,' ='
c      do i = 1, ncomps(nivel)
c      write(*, *) (jacx(i, j),j=1,n)
c      end do
c      write(*, *)
c---------------------------------------------------------------------

c      read(*, *) nada

c  Compute New y
      do i = 1, ncomps(nivel)
      ynew(i) = b(i, nivel)
      do j = 1, ncpan
      ynew(i) = ynew(i) + a(i, j, nivel)*y(j)
      end do
      end do

c  Compute y at nivel and derivative of the corresponding row

      do i = 1, ncomps(nivel)
      call redneupsi(npsi, ynew(i), y(i), dpsi, nivel, depth) 
      do j = 1, n
      jacx(i, j) = jacx(i, j)*dpsi
      end do

      do j = 1, ncolult
      jaca(i, j) = jaca(i, j)*dpsi
      end do

      end do



c      write(*, *)' Updated y at level ', nivel, ' = '
c      write(*, *)(y(i), i=1, ncomps(nivel))

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c      read(*, *) nada

c      write(*, *)

      end do

c This end do corresponds to "do nivel=1,depth"
c---------------------------------------------------------------------

c      write(*, *)' Jacobian jaca is complete. Columns:', ncolult,
c     *  ' Rows:', m

c      do i = 1, m
c      write(*, *)' Row ', i 
c      write(*, *)(jaca(i, j),j=1,ncolult)
c      end do



c      write(*, *)' m = ', m,' ncomps(depth)=', ncomps(depth),' iguales'
      do i = 1, m
      f(i) = y(i)
      end do

c      write(*, *)' F =', (f(i),i=1,m)

c      read(*, *) nada
c      write(*, *)' We return from redneu'
c      read(*, *) nada
      

      return
      end


      
      subroutine redneupsi(npsi, x, psi, dpsi, nivel, depth)
      implicit none
      double precision x, psi, dpsi, z
      integer npsi, nivel, depth

      if(npsi.eq.1) then
      psi = x
      dpsi = 1.d0
      return
      endif
      if(npsi.eq.2) then
c RELU
      if(x.le.0.d0) then
      psi = 0.d0
      dpsi = 0.d0
      else
      psi = x
      dpsi = 1.d0
      endif
      return
      endif

c RELU except for nivel = depth
      if(npsi.eq.3) then
      if(x.le.0.d0) then
      psi = 0.d0
      dpsi = 0.d0
      else
      psi = x
      dpsi = 1.d0
      endif
      if(nivel.eq.depth) then
      psi = x
      dpsi = 1.d0
      endif
      return
      endif


c RELU smoothed
      if(npsi.eq.4) then
      z = dsqrt(x*x + 1.d-3)
      psi = (z + x)/2.d0 
      
      dpsi =  ((x/z) + 1.d0)/2.d0
      return
      endif

c  NN-sin
      if(npsi.eq.5) then
      psi = dsin(x)
      dpsi = dcos(x)
      return
      endif

 
c  NN- 2sin
      if(npsi.eq.6) then
      psi = 2.d0 * dsin(x)
      dpsi = 2.d0* dcos(x)
      return
      endif        

      end 





                       
      subroutine rondo(seed, x)

c   Random between -1 and 1
c
C     This is the random number generator of Schrage:
C
C     L. Schrage, A more portable Fortran random number generator, ACM
C     Transactions on Mathematical Software 5 (1979), 132-138.

      double precision seed, x 

      double precision a,p,b15,b16,xhi,xalo,leftlo,fhi,k
      data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/

      xhi= seed/b16
      xhi= xhi - dmod(xhi,1.d0)
      xalo= (seed-xhi*b16)*a
      leftlo= xalo/b16
      leftlo= leftlo - dmod(leftlo,1.d0)
      fhi= xhi*a + leftlo
      k= fhi/b15
      k= k - dmod(k,1.d0)
      seed= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
      if (seed.lt.0) seed = seed + p
      x = seed*4.656612875d-10

      x = 2.d0*x - 1.d0

      return

      end 
            

      subroutine randin(seed, menor, mayor, sorteado)
      implicit none
      double precision seed, z
      integer menor, mayor, sorteado
c  Computes a random integer between menor and mayor
c  menor must be less than or equal to mayor
1     call rando(seed, z)
      z = dfloat(menor) + z*dfloat(mayor + 1 -menor)
      sorteado = z
      if(sorteado.lt.menor.or.sorteado.gt.mayor) go to 1
      return
      end 




                       
      subroutine rando(seed, x)

C     This is the random number generator of Schrage:
C
C     L. Schrage, A more portable Fortran random number generator, ACM
C     Transactions on Mathematical Software 5 (1979), 132-138.

      double precision seed, x 

      double precision a,p,b15,b16,xhi,xalo,leftlo,fhi,k
      data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/

      xhi= seed/b16
      xhi= xhi - dmod(xhi,1.d0)
      xalo= (seed-xhi*b16)*a
      leftlo= xalo/b16
      leftlo= leftlo - dmod(leftlo,1.d0)
      fhi= xhi*a + leftlo
      k= fhi/b15
      k= k - dmod(k,1.d0)
      seed= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
      if (seed.lt.0) seed = seed + p
      x = seed*4.656612875d-10

      return

      end 
             
 	subroutine chole (n, a, posdef, add, nlin)
	implicit double precision (a-h,o-z)
	logical posdef
	dimension a(nlin, n)
	posdef = .true.


	if(a(1,1) .le. 0.d0) then
	posdef = .false.
	return
	endif

	a(1,1) = dsqrt(a(1,1))
	if(n.eq.1)return
	do 1 i=2,n
	do 2 j=1,i-1
	z = 0.d0
	if(j.gt.1)then
	do 3 k=1,j-1
	z = z + a(i,k) * a(j,k)
3     continue
	endif
	a(i,j) = (a(i,j) - z)/a(j,j)
2	continue
	z = 0.d0
	do 4 j=1,i-1
	z = z + a(i,j)**2
4     continue
	if( a(i,i) - z .le. 0.d0) then
	posdef = .false.
	add = z - a(i,i)
	return
	endif

	a(i,i) = dsqrt( a(i,i) - z )
1	continue
	return
	end



 	subroutine sicho (n, a, x, b, aux, nlin)
	implicit double precision (a-h,o-z)
	dimension a(nlin,n),x(n),b(n),aux(n)
	aux(1) = b(1)/a(1,1)
	if(n.gt.1)then
	do 1 i=2,n
	z = 0.d0
	do 2 j=1,i-1
	z = z + a(i,j)*aux(j)
2     continue
	aux(i) = (b(i) - z) / a(i,i)
1	continue
	endif
	x(n) = aux(n)/a(n,n)
	if(n.eq.1)return
	do 3 i=n-1,1,-1
	z = 0.d0
	do 4 j=i+1,n
	z = z + a(j,i)*x(j)
4     continue
	x(i) = (aux(i) - z)/a(i,i)
3	continue
	return
	end


 



      
