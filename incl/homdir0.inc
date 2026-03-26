C
C Le FICHIER homdir0.inc est l'ancienne version de $MEFISTO/incl/homdir.inc
C Il contient le NOM DU REPERTOIRE SUPPORT DU LOGICIEL MEFISTO
      CHARACTER*18   HOMDIR
      PARAMETER     (HOMDIR='/usr/local/mefisto')
cccC
cccC     ATTENTION: Depuis le 11/11/2009 le fichier homdir.inc est ecrase
cccC     lors de l'execution de cbl_tout
cccC     dans instalsource.bash par le nom actuel
cccC     du repertoire d'installation du source de Mefisto
cccC     Par exemple:
cccC     /home/alain/mefistosource/incl/homdir.inc
ccc      CHARACTER*25 HOMDIR
ccc      PARAMETER (HOMDIR='/home/monpc/mefistosource')
