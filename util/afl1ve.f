      SUBROUTINE AFL1VE( TITRE, NBC, VECT )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: AFFICHER LE MINIMUM, MOYENNE, MAXIMUM des COMPOSANTES DE VECT(NBC)
C ---- NON NULLES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR   Janvier 2012
C MODIFS: ALAIN PERRONNET Saint Pierre du Perray                Mai 2023
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      CHARACTER*(*)    TITRE
      DOUBLE PRECISION VECT(NBC), S, VMX, VMI, V

      NMI = 1
      NMX = 1
      VMX   = VECT(1)
      VMI   = VECT(1)
      S     = 0D0
      NBCN0 = 0

      DO N=1,NBC

         V = VECT(N)

         IF( V .LT. VMI ) THEN
            VMI = V
            NMI = N
         ENDIF

         IF( V .GT. VMX ) THEN
            VMX = V
            NMX = N
         ENDIF

         IF( V .NE. 0D0 ) THEN
            NBCN0 = NBCN0 + 1
            S = S + V
         ENDIF

      ENDDO


      IF( LANGAG .EQ. 0 ) THEN

ccc         print *,TITRE,NBCN0,' coefficients NON NULS parmi ', NBC
ccc         print *,TITRE,' Min Vect(n)=', VMI,' ncMIN=',NMI
ccc         IF( NBCN0 .EQ. 0 ) NBCN0 = 1
ccc         print *,TITRE,' Moy Vect(n)=', S/NBCN0
ccc         print *,TITRE,' Max Vect(n)=', VMX,' ncMAX=',NMX

         IF( NBCN0 .EQ. 0 ) NBCN0 = 1
         print *,TITRE,' Min Vect(n)=', VMI,' Moy Vect(n)=', S/NBCN0,
     %                 ' Max Vect(n)=', VMX

      ELSE

ccc         print *,TITRE,NBCN0,' NON ZERO coefficients among ',NBC
ccc         print *,TITRE,' Min  Vect(n)=', VMI,' ncMIN=',NMI
ccc         IF( NBCN0 .EQ. 0 ) NBCN0 = 1
ccc         print *,TITRE,' Mean Vect(n)=', S/NBCN0
ccc         print *,TITRE,' Max  Vect(n)=', VMX,' ncMAX=',NMX

         IF( NBCN0 .EQ. 0 ) NBCN0 = 1
         print *,TITRE,' Min Vect(n)=', VMI,' Mean Vect(n)=', S/NBCN0,
     %                 ' Max Vect(n)=', VMX

      ENDIF

      RETURN
      END
