!
!	PROGRAMA DE CALCULO DO TENSOR DESMAGNETIZACAO DE
!	POLIGONOS IRREGULARES (SLANTED)
!
!	O arquivo de entrada (IN_demag3D.dat) deve conter:
!	px(1),py(1): coordenadas do canto superior esquerdo (nm)
!	px(2),py(2): coordenadas do canto superior direito  (nm)
!	px(3),py(3): coordenadas do canto inferior direito  (nm)
!	px(4),py(4): coordenadas do canto inferior esquerdo (nm)
!	t	   : espessura (nm)
!
!	*************** RESTRICAO ****************
!	-------> px(1) e px(4) DEVEM ser iguais
!	-------> px(2) e px(3) DEVEM ser iguais
!	******************************************
!
!	O arquivo de saida (OUT_demag3D.dat) contem:
!	Linha 1: d1, d2, d3
!	Linha 2: d4, d5, d6
!	Linha 3: d7, d8, d9
!
!	OBS.: OS FATORES DE DESMAGNETIZACAO SAO OBTIDOS FAZENDO
!
!	      D11=d1/(4*Pi*V), D12=d2/(4*Pi*V), D13=d3/(4*Pi*V)
!	      D21=d4/(4*Pi*V), D22=d5/(4*Pi*V), D23=d6/(4*Pi*V)
!	      D31=d7/(4*Pi*V), D32=d8/(4*Pi*V), D33=d9/(4*Pi*V)
!
!	      ONDE V e DADO EM nm^3
!
!	A energia de desmagnetizacao (eV) e calculada como
!
!	Ud = K/2*(
!	          (d1*COS(phi)^2 + d5*SIN(phi)^2 + (d2+d4)*SIN(phi)*COS(phi))*SIN(theta)^2 +
!		  ((d3+d7)*COS(phi) + (d6+d8)*SIN(phi))*SIN(theta)*COS(theta) +
!		  d9*COS(theta)^2
!		 )
!	ONDE
!
!	K = mu0*Ms*Ms*JtoeV*ten**(-27)/(four*Pi)
!	Ms = 800 kA/m
!	mu0 = 4*Pi*10^-7 H/m
!	JtoeV = 6.242_lg*ten**(18)
!
!
!
!
MODULE precision_def
	IMPLICIT NONE
	INTEGER, PARAMETER :: lg=selected_real_kind(15,397)
	INTEGER, PARAMETER :: nmc=1000000  !1E6 !define o numero de amostras no metodo MC
	REAL(lg),PARAMETER :: zero=0._lg, one=1._lg, half=0.5_lg
END MODULE precision_def


MODULE dimensions_particle
	USE precision_def
	IMPLICIT NONE
	REAL(lg)  :: w,t
END MODULE dimensions_particle


MODULE vertices
	USE precision_def
	IMPLICIT NONE
	REAL(lg)  :: px(4), py(4)
END MODULE vertices


MODULE ab
	USE precision_def
	IMPLICIT NONE
	REAL(lg)  :: au, bu
	REAL(lg)  :: ad, bd
END MODULE ab



PROGRAM demag3D2
	USE precision_def
	USE dimensions_particle
	USE vertices
	USE ab
	IMPLICIT NONE
	REAL(lg)	:: ak(6), bk(6), ck(6)
	REAL(lg)	:: d(9)
	REAL(lg)	:: teta2,alfa2,beta2
	REAL(lg)	:: teta5,alfa5,beta5

	INTERFACE

		SUBROUTINE MC_demag(int_a, int_b, int_c)
		USE precision_def
		USE dimensions_particle
		USE vertices
		USE ab
		IMPLICIT NONE
		REAL(lg),INTENT(OUT) :: int_a(6),int_b(6),int_c(6)
		END SUBROUTINE MC_demag

	END INTERFACE


!	PRINT*,''
!	PRINT*,'---> RUNNING....'
!	PRINT*,''

100	FORMAT(E13.4,E13.4,E13.4)
  integer, allocatable :: seed(:)
	INTEGER :: n

	CALL RANDOM_SEED(size=n)
	allocate(seed(n))
	WRITE(*,*) n
	WRITE(*,*) 'SEED'
	call random_seed(get=seed)
	write (*, *) seed

	CALL RANDOM_SEED(put = seed+500)


	OPEN(unit=9,  file='IN_demag3D3.dat',  status='unknown')
	OPEN(unit=10, file='OUT_demag3D3.dat', status='unknown')

	READ(9,*)px(1),py(1)
	READ(9,*)px(2),py(2)
	READ(9,*)px(3),py(3)
	READ(9,*)px(4),py(4)
	READ(9,*)t
	CLOSE(unit=9)

	w  = px(2)-px(1)
	au = (py(2)-py(1))/w
	bu = py(1) - au*px(1)
	ad = (py(3)-py(4))/w
	bd = py(3) - ad*px(3)

	WRITE(*,*) w, au, bu,ad,bd

	CALL MC_demag(ak ,bk, ck)

	PRINT*, 'AK'
	WRITE(*,*) ak
	PRINT*, 'BK'
	WRITE(*,*) bk
	PRINT*, 'CK'
	WRITE(*,*) ck

	teta2 =  ATAN(au)
	alfa2 = -SIN(teta2)
	beta2 =  COS(teta2)

	WRITE(*,*) teta2,alfa2,beta2

	teta5 =  ATAN(ad)
	alfa5 =  SIN(teta5)
	beta5 = -COS(teta5)

	WRITE(*,*) teta5,alfa5,beta5

	d(1) = ak(1) + alfa2*ak(2) - ak(4) + alfa5*ak(5)
	d(2) = beta2*ak(2) + beta5*ak(5)
	d(3) = ak(3) - ak(6)

	d(4) = bk(1) + alfa2*bk(2) - bk(4) + alfa5*bk(5)
	d(5) = beta2*bk(2) + beta5*bk(5)
	d(6) = bk(3) - bk(6)

	d(7) = ck(1) + alfa2*ck(2) - ck(4) + alfa5*ck(5)
	d(8) = beta2*ck(2) + beta5*ck(5)
	d(9) = ck(3) - ck(6)

	PRINT*, 'D:'
	WRITE(*,*) d(1)
	WRITE(*,*) d(2)
	WRITE(*,*) d(3)
	WRITE(*,*) d(4)
	WRITE(*,*) d(5)
	WRITE(*,*) d(6)
	WRITE(*,*) d(7)
	WRITE(*,*) d(8)
	WRITE(*,*) d(9)

	WRITE(10,100)d(1),d(2),d(3)
	WRITE(10,100)d(4),d(5),d(6)
	WRITE(10,100)d(7),d(8),d(9)

	CLOSE(unit=10)


END PROGRAM demag3D2



	SUBROUTINE MC_demag(int_a, int_b, int_c)
	USE precision_def
	USE dimensions_particle
	USE vertices
	USE ab
	IMPLICIT NONE
	REAL(lg),INTENT(OUT) :: int_a(6), int_b(6), int_c(6)
	REAL(lg) :: x(3), Dy
	REAL(lg) :: fa(6), fb(6), fc(6)
	REAL(lg) :: inta(6), intb(6), intc(6)
	REAL(lg) :: s(6)
	INTEGER  :: i0

	s(1)  = w*t*t*(py(2) - py(3))
	s(2)  = w*w*t*t
	s(3)  = w*w*t
	s(4)  = w*t*t*(py(1) - py(4))
	s(5)  = s(2)
	s(6)  = s(3)

	fa = zero
	fb = zero
	fc = zero


	DO i0=1,nmc

	   CALL xyz(x, Dy)
	   CALL fafbfc(inta, intb, intc)
	   fa = fa + inta
	   fb = fb + intb
	   fc = fc + intc

	ENDDO

  PRINT*, 'FA'
	WRITE(*,*) fa
	PRINT*, '\n\n FB'
	WRITE(*,*) fb
	PRINT*, '\n\n Fc'
	WRITE(*,*) fc
	PRINT*, '\n\n'

	int_a = s*fa/(REAL(nmc))
	int_b = s*fb/(REAL(nmc))
	int_c = s*fc/(REAL(nmc))



	CONTAINS


		FUNCTION frand(kmin, kmax)
		USE precision_def
		IMPLICIT NONE
		REAL(lg),INTENT(IN)  :: kmin,kmax
		REAL(lg)	     :: rnd0,frand
		CALL RANDOM_NUMBER(rnd0)
		!rnd0 = 0.50000000000000
		frand = kmin*(one - rnd0) + rnd0*kmax

		END FUNCTION frand


		SUBROUTINE fafbfc(int0a, int0b, int0c)
		USE precision_def
		IMPLICIT NONE
		REAL(lg),INTENT(OUT) :: int0a(6), int0b(6), int0c(6)
		REAL(lg)	     :: xp(3), int0(3), yminp, ymaxp
!		****S1****
		xp(1) = px(2)
		xp(2) = frand(py(3), py(2))
		xp(3) = frand(-half*t, half*t)
		CALL fu(xp,int0)
		int0a(1) = int0(1)
		int0b(1) = int0(2)
		int0c(1) = int0(3)
!		****S2****
		xp(1) = frand(px(1), px(2))
		xp(2) = au*xp(1)+bu
		xp(3) = frand(-half*t, half*t)
		CALL fu(xp,int0)
		int0a(2) = int0(1)
		int0b(2) = int0(2)
		int0c(2) = int0(3)
!		****S3****
		xp(1) = frand(px(1), px(2))
		yminp = ad*xp(1)+bd
		ymaxp = au*xp(1)+bu
		xp(2) = frand(yminp, ymaxp)
		xp(3) = half*t
		CALL fu(xp,int0)
		int0a(3) = int0(1)*(ymaxp - yminp)
		int0b(3) = int0(2)*(ymaxp - yminp)
		int0c(3) = int0(3)*(ymaxp - yminp)
!		****S4****
		xp(1) = px(1)
		xp(2) = frand(py(4), py(1))
		xp(3) = frand(-half*t, half*t)
		CALL fu(xp,int0)
		int0a(4) = int0(1)
		int0b(4) = int0(2)
		int0c(4) = int0(3)
!		****S5****
		xp(1) = frand(px(1), px(2))
		xp(2) = ad*xp(1)+bd
		xp(3) = frand(-half*t, half*t)
		CALL fu(xp,int0)
		int0a(5) = int0(1)
		int0b(5) = int0(2)
		int0c(5) = int0(3)
!		****S6****
		xp(1) = frand(px(1), px(2))
		yminp = ad*xp(1)+bd
		ymaxp = au*xp(1)+bu
		xp(2) = frand(yminp, ymaxp)
		xp(3) = -half*t
		CALL fu(xp,int0)
		int0a(6) = int0(1)*(ymaxp - yminp)
		int0b(6) = int0(2)*(ymaxp - yminp)
		int0c(6) = int0(3)*(ymaxp - yminp)
		END SUBROUTINE fafbfc


		SUBROUTINE xyz(x_0, Dy0)
		USE precision_def
		IMPLICIT NONE
		REAL(lg),INTENT(OUT) :: x_0(3), Dy0
		REAL(lg)             :: ymin0, ymax0
		x_0(1) = frand(px(1), px(2))
		ymin0 = ad*x_0(1)+bd
		ymax0 = au*x_0(1)+bu
		Dy0 = ymax0 - ymin0
		x_0(2) = frand(ymin0, ymax0)
		x_0(3) = frand(-half*t, half*t)
		END SUBROUTINE xyz


		SUBROUTINE fu(xp0, int_0)
		USE precision_def
		IMPLICIT NONE
		REAL(lg),INTENT(IN)  :: xp0(3)
		REAL(lg),INTENT(OUT) :: int_0(3)
		REAL(lg)	     :: den

		den = (&
					SQRT(&
					(x(1) - xp0(1))*(x(1) - xp0(1))+ &
			    (x(2) - xp0(2))*(x(2) - xp0(2))+ &
			    (x(3) - xp0(3))*(x(3) - xp0(3))&
					)&
					)**(-3)

		int_0(1) = Dy*(xp0(1) - x(1))*den
		int_0(2) = Dy*(xp0(2) - x(2))*den
		int_0(3) = Dy*(xp0(3) - x(3))*den
		END SUBROUTINE fu


	END SUBROUTINE MC_demag
