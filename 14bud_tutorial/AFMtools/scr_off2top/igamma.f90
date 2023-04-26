      program  main
        implicit none
        CHARACTER(100) :: num1char
       ! CHARACTER(100) :: num2char
        real*8 a,x,y,igamma
        CALL GET_COMMAND_ARGUMENT(1,num1char)
       ! CALL GET_COMMAND_ARGUMENT(2,num2char)
       ! read(num1char,*) a
        read(num1char,*) x      
        a=2.0d0/3.0d0
        print*, "Aequal",a        
!E = P(1)/R *
        !(1.d0-exp(-Rtmp)+P(2)**(1.d0/3.d0)*igamma(2.d0/3.d0, Rtmp))
        y=igamma(a, x)
        print*, y
      end
!
!
       Real*8 Function igamma(a, x) ! ¦£(a, x)=¦£(a) if x<=0 or a<=0
       !       ITMAX: the maximum allowed number of iterations
       !       EPS  : the relative accuracy
       !       FPMIN: a number near the smallest representable floating-point
       !       number
               integer:: n, ITMAX=100
               real*8 :: a, x, an, Tn, b, c, d, EPS=1.d3*epsilon(1.d0),FPMIN=1.d-30
               if(x<=0.d0 .or. a<=0.d0) then
                       igamma=gamma(a)
                       return
               end if
              
               igamma=0.d0
               if(x<a+1.d0) then ! use series representation
                       an=a; Tn=1.d0/a; igamma=Tn
                       do n=1, ITMAX
                               an = an+1.d0
                               Tn = Tn*x/an
                               igamma = igamma+Tn
                               if(abs(Tn)<abs(igamma)*EPS) then
                                       igamma = gamma(a)-igamma * exp(-x)*x**a
                                       return
                               end if
                       end do
               else ! use continued fraction representation
                       b=x+1.d0-a  ! modified Lentz's method with b0=0
                       c=1.d0/FPMIN
                       d=1.d0/b
                       igamma=d
                       do n=1, ITMAX
                               an=dble(-n)*(dble(n)-a)
                               b=b+2.d0
                               d=an*d+b; if(abs(d)<FPMIN) d=FPMIN
                               c=b+an/c; if(abs(c)<FPMIN) c=FPMIN
                               d=1.d0/d
                               Tn=d*c
                               igamma = igamma*Tn
                               if(abs(Tn-1.d0)<EPS) then
                                       igamma = igamma * exp(-x) * x**a
                                       return
                               end if
                       end do
               end if
       End
