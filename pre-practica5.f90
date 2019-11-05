! ------------------------------------ Pràctica 5 --------------------------------------- !
! Autor: Javier Rozalén Sarmiento
! Grup: B1B
! Data: 05/11/2019
!
! Funcionalitat: es programen els mètodes de "sampleig" de densitats de probabilitat de canvi
! de variable i d'acceptació i rebuig, i es posen a prova amb algunes funcions.
!
! Comentaris: tot i que el meu grup és el B1B, tal i com vaig quedar amb en Bruno, entrego
! aquesta pràctica el mateix dia 05/11/2019 però a la sessió de les 15h-17h.

program practica5
    implicit none
    double precision pi,M,a,b,fun
    double precision, allocatable :: xnums(:)
    integer ndat,i
    external fun
    common/cts/pi
    pi=acos(-1.d0)

    ! -------------------------------- Execici 1 --------------------------------------- !
    ndat=2000
    a=0.d0
    b=pi
    M=1.d0
    allocate(xnums(ndat))
    call acceptrebuig(ndat,xnums,a,b,M,fun)
    do i=1,40
        print*,xnums(i)
    enddo
end program practica5

subroutine acceptrebuig(ndat,xnums,a,b,M,funci)
    implicit none
    double precision xnums(ndat),a,b,M,funci,x1,x2,p,x
    integer ndat,iseed,counter
    iseed=20034276  
    counter=0
    call srand(iseed)
    ! Algorisme
    1   x1=rand()
    x2=rand()
    x=(b-a)*x1+a
    p=M*x2
    if (funci(x).ge.p) then 
        xnums(counter)=x
        counter=counter+1
    else
        goto 1
    endif
    return
end subroutine acceptrebuig

double precision function fun(x)
    implicit none
    double precision pi,x
    common/cts/pi
    fun=(12.d0/(pi*(2.d0*pi**2.d0-3.d0)))*x**2.d0*(sin(x))**2.d0
    return
end function fun

SUBROUTINE HISTOGRAMA(NDAT,XDATA,XA,XB,NBOX,XHIS,VHIS,ERRHIS,BOXSIZE,IERR)
    IMPLICIT NONE
    ! INPUT/OUTPUT VARIABLES
    INTEGER NDAT,NBOX
    DOUBLE PRECISION XDATA(NDAT),XA,XB
    DOUBLE PRECISION XHIS(NBOX),VHIS(NBOX),ERRHIS(NBOX)
    INTEGER IERR

    INTEGER I,IBOX,ICOUNT
    DOUBLE PRECISION BOXSIZE

    IF (XA.GE.XB) THEN 
        IERR=1
        RETURN
    ENDIF
    ! BOX SIZE
    BOXSIZE=(XB-XA)/NBOX

    ! COUNTS NUMBER OF POINTS WITHIN THE INTERVAL XA,XB
    ICOUNT=0

    ! SETS ALL TO ZERO
    DO I=1,NBOX
        VHIS(I)=0
        ERRHIS(I)=0
    ENDDO

    ! WE RUN THROUGH THE DATASET
    DO I=1,NDAT
    ! CHECKS IF DATA LIES WITHIN XA,XB
    IF (XDATA(I).GE.XA.AND.XDATA(I).LE.XB) THEN 
        IBOX=INT((XDATA(I)-XA)/BOXSIZE)+1
    ! PUTS XB INTO THE LAST BOX, IF NEEDED
        IF (IBOX.EQ.NBOX+1) IBOX=NBOX 

            VHIS(IBOX)=VHIS(IBOX)+1
            ICOUNT=ICOUNT+1
    ENDIF
    ENDDO

    IF (ICOUNT.EQ.0) THEN 
       IERR=2
       RETURN
    ENDIF

    IERR=0
    PRINT*,"ACCEPTED:",ICOUNT," OUT OF:",NDAT

    DO I=1,NBOX
    ! CENTRAL VALUE OF THE BAR
        XHIS(I)=XA+BOXSIZE/2.D0+(I-1)*BOXSIZE
    !  ERROBAR, STANDARD DEVIATION OF CORRESPONDING BINOMIAL
        ERRHIS(I)=SQRT(VHIS(I)/ICOUNT*(1.D0-VHIS(I)/ICOUNT))/BOXSIZE / SQRT(DBLE(ICOUNT))
    ! NORMALIZED VALUE OF THE BAR
    VHIS(I)=VHIS(I)/ICOUNT/BOXSIZE
    ENDDO
END SUBROUTINE HISTOGRAMA

