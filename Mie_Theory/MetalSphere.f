c
c SKG program to create a script file for running the multilayer code
c MODIFIED by JJF

c Five Layers:  
c (1) Au core 
c (2) Au reduced conductivity layer
c (3) OAm Surfactant layer
c (4) Fe Layer
c (5) Fe3O4 Layer
c (6) Fe2O3 Layer

c radius  = radius of gold core - dR
c dR      = 1 atomic layer
c dFe     = thickness of Fe middle shell
c dM      = thickness of Fe3O4 second shell
c dFeO    = Thickness of Fe2O3 outter shell 
c r = radius + dFe + dM + dFeO

       implicit real*8(a-h,o-z)
       complex*16 eps(6),m(6),ii,perm,RI
       real*8 mreal(6),mimag(6),r(6),x(6)
c RI for various materials/solvents.
       ii=(0.d0,1.d0)
       water=1.33d0
       air=1.0d0
       hexane=1.37d0
       amine=1.46d0
       toluene=1.5
       glass=1.6
       tio2=2.5
c specify these values (nm or no units), run code and generate script:
c edit as appropriate.
c
c ********************************************************
c ********************************************************
c diameter of Au core:
        read(5,*)diameter

c radius of shell 1 (QD)
        read(5,*)refmed

        radius=diameter/2.d0
c thickness of lower conductivity layer:
c n.b. one layer of silver might be 2.5 angstroms or 0.25 nm thick ...
c RI of surrounding Ag layer

c ********************************************************
c ********************************************************
        write(11,*)refmed
        ii=(0.d0,1.d0)
        pi=4.d0*datan(1.d0)
        write(6,*)'#! /bin/sh'
c
c input lines to scattnlay are 
c -l NumberLayers  x1  real(m1)  imag(m1)  x2  real(m2) imag(m2) .... 
c Here, xi = 2 pi Nm ri /lambda,   Nm = surrounding medium RI, ri = radius of layer and
c mi = Ni/Nm.
c
        nlayer=1
c
        r(1)=radius
c
c       rlambda will start at 200 nm
        rdlambda=800.d0/1000.d0

        open(unit=8,file=
     >  'DIEL/Pt_Palik.txt')

        do jj=1,1778
        read(8,*)rlambda,epsre,epsim
        eps(1)=epsre+ii*epsim
        m(1)=cdsqrt(eps(1))/refmed


        do i=1,nlayer
        x(i)=2.d0*pi*refmed*r(i)/rlambda
        mreal(i)=dreal(m(i))
        mimag(i)=dimag(m(i))
        enddo

        write(6,1066)nlayer,
     >  x(1),mreal(1),mimag(1),rlambda
c
c
        enddo
c
        close(8)
 1066   format(1x,'SCAT_CODE/scattnlay -l',i4,3e12.5,' -c ',e12.5)
c
c
        end
