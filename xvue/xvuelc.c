/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* BUT : INTERFACE ENTRE FORTRAN ET X11 VIA LE LANGAGE  POUR FAIRE                */
/* ----- DES TRACES GRAPHIQUES SOUS X11 PORTABLES                                 */
/*       INTERFACE POUR HABILLER Mefisto            SANS MOTIF ni TOOLKIT         */
/*       IMPRESSION DES TRACES EN POSTSCRIPT                                      */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* AUTEUR : PERRONNET ALAIN    ANALYSE NUMERIQUE UPMC  PARIS       NOVEMBRE 1994  */
/* AUTEUR : DOURSAT CHRISTOPHE ANALYSE NUMERIQUE UPMC  PARIS       NOVEMBRE 1994  */
/* MODIFS : PERRONNET ALAIN LJLL UPMC PARIS & St PIERRE DU PERRAY  OCTOBRE  2013  */
/*................................................................................*/  
#include <limits.h>    /*limites min max int long real ... */
#include <unistd.h>
#include <stdio.h>	
#include <string.h>		
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <setjmp.h>
#include <ctype.h>   /* Ajout 20/07/2012 pour toupper  */

#include <sys/types.h>
#include <sys/time.h>

/*
   Pour retrouver les fichiers X11, charger le paquet libx11-dev
*/

/*
#include <X11/Xlib.h>
#include <X11/Xatom.h>
#include <X11/Xutil.h>
*/

/*
   Les 3 lignes suivantes sont necessaires pour UBUNTU 10.10 ou Mint
   car non retrouvees dans <X11/...>  (pas necessaires pour 11.04)
   Attention: taper gcc avec -I. pour les retrouver!
              gcc -O -I. -c xvue/xvuelc.c -o xvue/xvuelc.o
*/

#include "./incl/Xlib.h"
#include "./incl/Xatom.h"
#include "./incl/Xutil.h"

/* fin 2007 include X11 devenus inutiles
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Shell.h>
*/

/* =============================================================================*/
/*    PARTIE DEPENDANTE DES APPELS FORTRAN DE FONCTIONS C                       */
/*        A AJUSTER AU MIEUX SELON LA VERSION LOCALE                            */ 
/*                                                                              */
/*    macro  proc( )    Procedure necessaire pour generer un _                  */
/*                     derriere le nom des fonctions C afin                     */
/*                     qu'elles puissent etre appelees en FORTRAN               */
/*                  => GRANDE DEPENDANCE AVEC C et FORTRAN et LINK local        */
/* =============================================================================*/ 
#ifdef __GNUC__
#define  proc(x) x##_
/*     a choisir sur PC_LINUX DEC SUN SGI et HP */
#else
#define  proc(x) x/**/_
/*  parfois proc(x) x    SUFFIT !  (cas IBM) */
#endif

/* LE CHOIX pour cette STATION  */
#undef  proc
#define proc(x) x##_

/* =============================================================================*/

#define min(A, B) ((A) < (B) ? (A) : (B))       /* MIN et MAX de 2 valeurs */
#define max(A, B) ((A) > (B) ? (A) : (B))

#define MAX_SHORT 0xFFFF   /* Plus grande valeur d'un entier short */

/* =============================================================================*/
/*         Les constantes nommees pour X et la fenetre mefisto                  */
/* =============================================================================*/
#define X_COIN_WINDOW 0  /* POSITION du COIN SUPERIEUR GAUCHE DE LA FENETRE X   */
#define Y_COIN_WINDOW 0

#define	CMAPSIZE     256  /* Nombre de couleurs mefisto pour 8 Plans 256=2**8 */

/* =============================================================================*/
/* La ZONE DES DONNEES pour X et le graphisme                                   */
/* =============================================================================*/
static  int       langage;       /* 0=>Francais, 1=>Anglais */

static  Display  *display_mef; /* la structure display est en fait le SERVEUR X */

static  GC gc_mef;             /* contexte graphique pour la fenetre mef*/

static  Window fenetre_mef0, fenetre_mef; /*caracteristiques de la fenetre actuelle*/

static  XSetWindowAttributes  attributes;

static  float   red[CMAPSIZE], 
              green[CMAPSIZE], 
               blue[CMAPSIZE];  /* rouge vert bleu de la palette des couleurs */

static unsigned long norgb[CMAPSIZE];
              /* numero dans la palette physique de la couleur i dans la palette
                 des couleurs de mefisto
              en PseudoColor: c'est l'identite norgb[i]=i
                              la i_eme couleur de la palette mefisto est celle 
                              de la palette physique
              en TrueColor  : c'est le no dans la palette physique en lecture seule
                              norgb[i]=fonction de (rouge, vert, bleu)
                              cette valeur est calculee par XAllocColor
                              qui determine a partir de RGB la couleur la plus 
                              proche parmi les couleurs de la palette physique */

unsigned long foreground, background; /* couleur de trace et couleur du fond */

static	int       screen_mef;    /* caracteristiques de l'ecran graphique */

static  Colormap  color_map0, color_map; 
                             /* palette initiale (En lecture SEULE)
                             et palette courante des couleurs (en LECTURE/ECRITURE) */

static Visual   *visual;  /* Choix du type d'utilisation de l'ecran PseudoColor,... */
                                                 
static XFontStruct *struc_police ;   /* structure de la police courante */
static int ascent_pol, descent_pol;  /* nombre de pixels au dessus et dessous de la ligne de base */
static XCharStruct mesure ;          /* qualifiants des dimensions d'une chaine de caracteres */
static int         lapxfe, lhpxfe;   /* nombre de pixels de la fenetre en largeur et hauteur */

int    nbpolices;    /* nombre de polices ou fontes de caracteres disponibles */        
char **listfonts;    /* pointe sur une liste de noms de fontes */
Font   font_id;      /* nom de la fonte actuelle */  
static long nbs1970; /* sert a dissocier les decimales des secondes depuis 1970 */
                     /* et surtout pour que 2 tableaux n'aient pas la meme date */

static Pixmap      mempx;    /* tous les traces y sont faits avant copie dans fenetre_mef */
static Pixmap      mempxsauvfen;  /* sauvegarde du contenu de la fenetre Mefisto
                                     avant affichage du rectangle du menu */

static  Pixmap  mempxaccro; /* definition du pixmap pour l'accrochage des poignees */
#define lmempxaccro 13      /* largeur en pixels du mempxaccro d'accrochage */
#define hmempxaccro 13      /* hauteur en pixels du mempxaccro d'accrochage */

/*=============================================================================*/ 
/*                ZONE DES DONNEES POSTSCRIPT                                  */
/*=============================================================================*/
/* parametres contenus dans les communs decrivant les tailles des divers menus */
/* dans gsmenu.inc  */
#define MXRECT  8   /* nombre maximal de rectangles menus                      */
#define NBCAME 80   /* nombre de caracteres d'une ligne du menu                */
#define MXLGME 48   /* nombre maximal de lignes du menu                        */
#define NBCAHI 48   /* nombre de caracteres d'une ligne de l'historique        */
#define MXLGHI 16   /* nombre maximal de lignes de l'historique                */
#define NBCADO 80   /* nombre de caracteres d'une ligne de la documentation    */
#define MXLGDO 31   /* nombre maximal de lignes de la documentation            */
#define NBCAIN 48   /* nombre de caracteres d'une ligne de l'invite            */
#define MXLGIN  2   /* nombre maximal de lignes de l'invite                    */
#define NBCALG 96   /* nombre de caracteres d'une ligne de saisie              */
#define MXLGSA  1   /* nombre maximal de lignes de la ligne de saisie          */
#define NBCAER 96   /* nombre de caracteres d'une ligne des erreurs            */
#define MXLGER 32   /* nombre maximal de lignes des erreurs                    */ 

/* dans lu.inc      */
#define NBCALI 96   /* nombre de caracteres par ligne lue                      */
#define MXKLG 256   /* nombre maximal de lignes de lignes lues                 */

#define NBCALIP10 106 /* NBCALI+10 pour ajout dans kliglue du texte Ligne lue: */
                                                                                
static int   lasopsc; /* sortie ou non du trace en postscript */
static int   modepsc ;
static FILE  *fpi, *fpo ; 
static float counb , courgb[3] ;
static int   palcourc ;
static int   xpixels , ypixels ;
static char* chaine[MXRECT] ;
static int   longchaine[MXRECT] ;
static char  format[255] ;
static int   menu ; 
static int   icolorm ; 
static int   n1core ;
static int   n1coel ;
static int   ndcoel ;
static int   ndcore ;
static int   n1coul ;
static int   ndcoul ;
static int   nbrcon, xinic, yinic, xcouc, ycouc ;
static int   iTe, iFa, ity, iep, iPo, ire , iRe, iel, iEl, iFP ;
static char  buf[512] , concat[512] , fontcour[512] ;
/******************************************************************************/

/* =============================================================================*/
/* La ZONE DES DONNEES pour X                                                   */
/* =============================================================================*/
#define NBLIGNESMENU 25 /* Nombre de lignes max avant le scroll sur le menu */

/* context pour les longjmp : un pour les menus ,un pour les warnings */
jmp_buf context_menu,context_error;

/* chaine_de_saisie : chaine de caracteres en sortie du C renvoyee vers FORTRAN
   taille_chaine de saisie : taille de la chaine effectivement saisie */
char *chaine_de_saisie;
int *taille_de_la_chaine_saisie;

/* kmenu : tableau de lignes de chaines de caracteres contenant le menu a traiter */
char kmenu[MXLGME][NBCAME];
 
/* nom en C du fichier a ouvrir pour la documentation */
char nom_du_fichier_de_documentation[72];

/* nom en C du home directory de mefisto */
char nom_homdir[72];

/* kerr : tableau contenant le message d'erreur devant etre affiche */
char kerr[NBCAER*MXLGER];

/* kinvite : tableau contenant le message d'invite devant etre affiche */
char kinvite[NBCAIN*MXLGIN];

/* kliglue : tableau contenant le message des lignes lues devant etre affiche */
char kliglue[NBCALIP10*MXKLG];

/*============================================================================*/
/*             LES FONCTIONS ET PROCEDURES C                                  */
/*============================================================================*/

void proc(languemefisto)( int * langue )
{
  if ((fpo = fopen("/usr/local/mefisto/td/m/anglais","r"))!=NULL)
  {
      *langue=1;    /* langue anglaise */
      fclose(fpo) ;
      fpo = NULL ; 
  }
  else
  {
    *langue=0;      /* langue francaise */
  }
}


void * proc(dctnmc) ( int *nboctets )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT : Allouer dynamiquement un tableau de nboctets
 ----- Recuperer son adresse memoire dans mc et 0 si pas d'allocation
       sous la forme Tableau = dctnmc( NbOctets )
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AUTEUR : PERRONNET ALAIN  UPMC Laboratoire JLL  Paris             SEPTEMBRE 2005
12345X7..............................................................012345678*/
{
  return malloc( (size_t) (*nboctets) );
}

void  proc(dstnmc) ( void * mcoctets )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT : Desallouer un tableau allouer dynamiquement
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AUTEUR : PERRONNET ALAIN  UPMC Laboratoire JLL  Paris             SEPTEMBRE 2005
12345X7..............................................................012345678*/
{
  free( mcoctets );
  mcoctets = NULL;
  return;
}

void  proc(nomrepmefisto) ( char *chaine, int *size )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT : nomrepmefisto : passe le nom du repertoire de mefisto du FORTRAN au C
 -----
 ENTREES :
 ---------
 chaine : chaine a transformer 
 size   : taille de la chaine  
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AUTEUR : PERRONNET ALAIN     ANALYSE NUMERIQUE UPMC  PARIS          JANVIER 1995
12345X7..............................................................012345678*/ 
{
  int i;

  for ( i = 0; i < *size; i++ )
    nom_homdir[i] = chaine[i];
  nom_homdir[*size] = '\0';

}

void proc(xvinitgraphique) (void)  
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :   INITIALISER le graphisme en X (display_mef et screen_mef)
 -----   
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
 
{ 
  if ( (display_mef = XOpenDisplay((char *)NULL)) == (Display *)NULL)
  {
    printf("Initialisation Graphique X11 IMPOSSIBLE \n");
    printf("X11 CAN NOT BE INITIATED");
    exit(1);
  }
  screen_mef = DefaultScreen( display_mef ); /* L'ecran graphique */
  nbs1970 = 0;/* initialisation du nombre d'appels a la fonction secondes1970 */
}


void   proc(xtinit) ()
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT : xtinit : initialise le display et l'ecran
 -----
 remarque : cette routine doit etre appelle juste avant xvinit
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  proc(xvinitgraphique)();
}


void proc(xvpxecran)(int *xp, int *yp)   
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :   RECUPERER le nombre de pixels de la largeur et hauteur de l'ecran total
 -----   

 SORTIES :
 ---------
 xp      : nombre de pixels de la largeur de l'ecran
 yp      : nombre de pixels de la hauteur de l'ecran
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  *yp  = XDisplayHeight ( display_mef, screen_mef );
  *xp  = XDisplayWidth  ( display_mef, screen_mef );
}              


void proc(xvmmecran)(int *xmm, int *ymm)  
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :   RECUPERER la largeur et hauteur de l'ecran total en MILLIMETRES
 -----  

 SORTIES :
 ---------
 xmm     : nombre de MM de la largeur de l'ecran
 ymm     : nombre de MM de la hauteur de l'ecran
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  *ymm = XDisplayHeightMM ( display_mef, screen_mef );
  *xmm = XDisplayWidthMM  ( display_mef, screen_mef );
  //printf("xvmmecran avant: ecran xmm=%d ymm= %d\n", *xmm, *ymm);
  //*ymm = 180;  //sur portable vaio ap
  //*xmm = 287;  //sur portable vaio ap
  printf("xvmmecran: screen xmm=%d ymm= %d\n", *xmm, *ymm);
 }  

void proc(xvCouleursImposees)( int n1coel, int *ndcoel,
                               float red[], float green[], float blue[] )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
 BUT :     Definir les couleurs imposees dans toute palette
 -----     ( procedure interne a xvue )   

 ENTREES :
 ---------
 n1coel  : numero de la premiere couleur imposee

 MODIFIES:
 ---------
 red, green, blue : les tableaux ROUGE VERT BLEU de la palette

 SORTIES :
 ---------
 ndcoel  : numero de la derniere couleur imposee
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                 AOUT 1997
12345X7..............................................................012345678*/ 
{
                 /* noir    */
  red[n1coel]     = 0.0 ;
  green[n1coel]   = 0.0 ;
  blue[n1coel]    = 0.0 ; 

                 /* rouge   */       
  red[n1coel+1]   = 1.0 ;
  green[n1coel+1] = 0.0 ;
  blue[n1coel+1]  = 0.0 ;

                 /* vert sombre */
  red[n1coel+2]   =  50./256. ;
  green[n1coel+2] = 200./256. ;
  blue[n1coel+2]  =  50./256. ;

                 /* bleu    */
  red[n1coel+3]   = 0.0 ;
  green[n1coel+3] = 0.0 ;
  blue[n1coel+3]  = 1.0 ;

                 /* cyan    */
  red[n1coel+4]   = 0.0 ;
  green[n1coel+4] = 0.8 ;
  blue[n1coel+4]  = 1.0 ;

                 /* jaune   */
  red[n1coel+5]   = 1.0 ;
  green[n1coel+5] = 1.0 ;
  blue[n1coel+5]  = 0.0 ;

                 /* magenta */
  red[n1coel+6]   = 1.0 ;
  green[n1coel+6] = 0.0 ;
  blue[n1coel+6]  = 1.0 ;

                 /* blanc   */
  red[n1coel+7]   = 1.0 ;
  green[n1coel+7] = 1.0 ;
  blue[n1coel+7]  = 1.0 ;

                 /* gris1 sombre bleute */
  red[n1coel+8]   =  80./256. ;
  green[n1coel+8] =  80./256. ;
  blue[n1coel+8]  = 100./256. ;

                 /* gris2 moyen bleute */
  red[n1coel+9]   = 150./256. ;
  green[n1coel+9] = 150./256. ;
  blue[n1coel+9]  = 178./256. ;

                 /* gris3 clair bleute */
  red[n1coel+10]   = 220./256. ;
  green[n1coel+10] = 220./256. ;
  blue[n1coel+10]  = 256./256. ;

                 /* peachpuff  remplace  beige */
  red[n1coel+11]   = 256./256. ;  /* 250./256. ; */
  green[n1coel+11] = 218./256. ;  /* 240./256. ; */
  blue[n1coel+11]  = 185./256. ;  /* 216./256. ; */

                 /* orange   */
  red[n1coel+12]   =  1.0;
  green[n1coel+12] =  0.5;
  blue[n1coel+12]  =  0. ;

                 /* saumon */
  red[n1coel+13]   = 250./256. ;
  green[n1coel+13] = 128./256. ;
  blue[n1coel+13]  = 114./256. ;

                 /* rose */
  red[n1coel+14]   = 1.0 ;
  green[n1coel+14] = 190./256. ;
  blue[n1coel+14]  = 206./256. ;

                 /* turquoise */
  red[n1coel+15]   =  74./256. ;
  green[n1coel+15] = 250./256. ;
  blue[n1coel+15]  = 160./256. ;

  /* 16 COULEURS IMPOSEES */
  *ndcoel = n1coel + (int) 15 ;
}

void proc(xvColormapToRGB)( Colormap color_map, 
                            float r[], float g[], float b[], int nbcolor)    
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :   TRANSFERER les nbcolor couleurs de la palette color_map dans les
 -----   tableaux r, v, b de C  dans la classe PseudoColor
        ( procedure interne a xvue )   

 ENTREE :
 --------
 nbcolor : nombre de couleurs dans la palette des couleurs a charger

 SORTIES :
 ---------
 r, g, b : les tableaux ROUGE VERT BLEU de la palette en C
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                 AOUT 1997
12345X7..............................................................012345678*/ 
{
  /* XColor          rgb_def, hardware_def;
     unsigned short  maxvalue; */

  XColor    colorcell_defs;
  int       i;

  /* XLookupColor (display_mef, color_map, "white", &rgb_def, &hardware_def);
  maxvalue = hardware_def.red;      et maxvalue remplace les MAX_SHORT ci-dessous */

  for ( i = 0 ; i < nbcolor ; i++ )
  {
    colorcell_defs.pixel = (unsigned long) i;
    /* Recuperation du rouge vert bleu de la couleur i de la colormap */
    XQueryColor (display_mef, color_map, &colorcell_defs);

    r[i] = ( (float) colorcell_defs.red   / (float) MAX_SHORT );
    g[i] = ( (float) colorcell_defs.green / (float) MAX_SHORT );
    b[i] = ( (float) colorcell_defs.blue  / (float) MAX_SHORT );
  }
}     


void proc(xvStockeRGBtoColormap)( int nbcells,
                                  float red[], float green[], float blue[],
                                  unsigned long norgb[], Colormap color_map )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
 BUT :     IMPOSER dans la palette des couleurs color_map modifiable
 -----     les couleurs red green blue tableaux C
           ( procedure interne a xvue )   

 ENTREES :
 ---------
 nbcells : nombre de couleurs dans la palette des couleurs a charger
 red, green, blue : les tableaux ROUGE VERT BLEU de la palette en C ou en FORTRAN

 SORTIES :
 ---------
 norgb   : numero des couleurs de la palette mefisto dans la palette physique
 color_map : recoit les nbcells couleurs
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                 AOUT 1997
12345X7..............................................................012345678*/ 
{
  unsigned long  index;
  XColor         colorcell_defs[CMAPSIZE];

  if ( (visual->class != TrueColor) && (visual->class != GrayScale) )  /* PseudoColor par exemple */
  {
    for ( index = 0 ; index < nbcells ; index++ )
    {
    norgb[index] = (unsigned long) index; /* valeur ici imposee */
    colorcell_defs[index].pixel = (unsigned long) index; 
    colorcell_defs[index].red   = (unsigned short) ( red[index]   * (float) MAX_SHORT );
    colorcell_defs[index].green = (unsigned short) ( green[index] * (float) MAX_SHORT );
    colorcell_defs[index].blue  = (unsigned short) ( blue[index]  * (float) MAX_SHORT );
    colorcell_defs[index].flags = DoRed | DoGreen | DoBlue; 
    colorcell_defs[index].pad   = '.';
    }  
    /* Stockage des couleurs mefisto dans la colormap lecture/ecriture NON PARTAGEABLE */
    XStoreColors ( display_mef, color_map, colorcell_defs, nbcells );
  }
  else   /* TrueColor ou GrayScale */
  {
    for ( index = 0 ; index < nbcells ; index++ )
    {
    /* colorcell_defs[index].pixel = (unsigned long) index; */
    colorcell_defs[index].red   = (unsigned short) ( red[index]   * (float) MAX_SHORT );
    colorcell_defs[index].green = (unsigned short) ( green[index] * (float) MAX_SHORT );
    colorcell_defs[index].blue  = (unsigned short) ( blue[index]  * (float) MAX_SHORT );
    colorcell_defs[index].flags = DoRed | DoGreen | DoBlue; 
    colorcell_defs[index].pad   = '.';

    /* Initialisation de colorcell_defs[index].pixel par recherche X de la couleur
       RGB la plus proche dans la palette physique  palette lecture seule PARTAGEABLE */
    XAllocColor ( display_mef, color_map, &colorcell_defs[index] );
    norgb[index] = colorcell_defs[index].pixel;  /* valeur ici retrouvee */
    }
  }
}

void proc(initaccrochage) (void)
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :   Initialise le mempxaccro contenant le carre de visualisation d'un accrochage
 -----
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AUTEUR : PERRONNET ALAIN    ANALYSE NUMERIQUE UPMC  PARIS           OCTOBRE 2000
12345X7..............................................................012345678*/  
{  
  XGCValues gcvalues;
  unsigned long  blanc, noir, colour;

  /* indirection necessaire pour traiter a la fois PseudoColor et TrueColor */
  /* colour est le numero de la couleur dans mefisto dans la palette physique */ 
  /*printf("dans initaccrochage\n");
    for (k=0; k<256; k++) printf("couleur %d = %d\n",k,norgb[k]);*/

  blanc = WhitePixel ( display_mef, screen_mef );
  noir  = BlackPixel ( display_mef, screen_mef );

  /* Creation du pixmap d'un carre pour copie dans la fenetre de Mefisto*/
  mempxaccro = XCreatePixmap( display_mef, fenetre_mef, lmempxaccro, hmempxaccro,
	   	           DefaultDepth( display_mef, screen_mef ) );

  /*colour=norgb[background];*/
  colour=blanc;
  /*printf("dans initaccrochage background= %d colour= %d\n", background, colour);*/
  XSetForeground( display_mef, gc_mef, colour ); /* pixmap est mis a blanc */
  XFillRectangle( display_mef, mempxaccro, gc_mef, 0, 0, lmempxaccro, hmempxaccro );

  /* 3 epaisseurs de trait */
  gcvalues.line_width = 3;
  XChangeGC( display_mef, gc_mef, GCLineWidth, &gcvalues );

  /* cercle */
  /* colour = norgb[1];
  printf( "dans initaccrochage cercle rouge=1 colour= %d\n", colour );
  XSetForeground (display_mef, gc_mef, colour);
  XDrawArc( display_mef, mempxaccro, gc_mef, 2, 2, lmempxaccro-3, hmempxaccro-3, 0, 360*64); */
  
  /* dessin dans le pixmap du carre */
  /*colour=norgb[foreground];*/
  colour=noir;
  /*printf("dans initaccrochage foreground= %d colour= %d\n", foreground, colour);*/
  XSetForeground( display_mef, gc_mef, colour );
  XDrawRectangle( display_mef, mempxaccro, gc_mef, 1, 1, lmempxaccro-3, hmempxaccro-3 );

  /* retour a 1 epaisseur de trait */
  gcvalues.line_width = 1;
  XChangeGC(display_mef, gc_mef, GCLineWidth, &gcvalues);
 }

void proc(xvinfo)( int *ix, int *iy, int *maxfonts,
              int *n1coref, int *ndcoref, int *n1coelf,
              int *ndcoelf, int *n1coulf, int *ndcoulf, int *nbcolo,
              char namefonts[][256], int nbchar[], int *nbfonts, int *visuclass )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :   DEFINIR la largeur et hauteur PIXELS de la FENETRE A OUVRIR
 -----   a partir du coin superieur gauche de l'ecran 
         DEFINIR la palette des couleurs
         Recuperer en FORTRAN les caracteristiques des fontes de caracteres
          
 ENTREES :
 --------- 
 ix       : nombre de pixels de la largeur de la fenetre a ouvrir
 iy       : nombre de pixels de la hauteur de la fenetre a ouvrir 
 maxfonts : nombre maximal de fontes de caracteres stockables

 SORTIES :
 ---------   
 n1coref   : numero dans la palette des couleurs de la premiere couleur reservee
 ndcoref   : numero dans la palette des couleurs de la derniere couleur reservee    
 n1coelf   : numero dans la palette des couleurs de la premiere couleur elementaire
 ndcoelf   : numero dans la palette des couleurs de la derniere couleur elementaire
 n1coulf   : numero dans la palette des couleurs de la premiere couleur modifiable
 ndcoulf   : numero dans la palette des couleurs de la derniere couleur modifiable
 nbcolo    : nombre de couleurs disponibles de la palette des couleurs
 namefonts : noms des fontes disponibles 
 nbchar    : nombre de caracteres du nom de chacune des fontes disponibles  
 nbfonts   : nombre -1 de fontes disponibles (No de fonte de 0 a nbfonts)
 visuclass : 0 PseudoColor
             1 GrayScale
             2 DirectColor
             3 StaticColor
             4 StaticGray
             5 TrueColor
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
    int     x, y, i, k;
    /*char   *j;*/ 
    char    data[8];   /* Modif:  data[8] au lieu de data[7] 10/1/2010 */
 
    unsigned int  width, height, border_width;
      
    /*XFontStruct  *info;*/
    /*Font          font;*/

    /*XEvent        event;*/

    XGCValues     gcvalues;

    /* le no de la langue de Mefisto 0=>Francais, 1=>Anglais */
    proc(languemefisto)( &langage );

    x = (int) X_COIN_WINDOW;
    y = (int) Y_COIN_WINDOW;

    width        = (unsigned int) *ix;
    height       = (unsigned int) *iy;
    lapxfe       = width;
    lhpxfe       = height;
    border_width = (unsigned int) 5  ;

    fenetre_mef0   = RootWindow(    display_mef, screen_mef ); 
    visual    = DefaultVisual( display_mef, screen_mef );  

    relance:
    switch (visual->class)
    {
      /*-----------------------------------------*/
      /* Les 3 CLASSES avec PALETTES MODIFIABLES */
      /*-----------------------------------------*/

      case PseudoColor:  /* Classe des ecrans avec 8 plans ecrans couleurs => 256 couleurs  */
      /*---------------     Les palettes colormap sont en Lecture/ecriture non partageables */
      /*                    Le no couleur mefisto = celui de la palette physique            */
      /*                    Tout pixel de l'ecran a un no de couleur dans la palette        */
      /*                    auquel correspond dans la palette 3 valeurs rouge vert bleu tracees */
      {
	if( langage )
	  printf(" Class PseudoColor\n");
	else
	  printf(" Classe PseudoColor\n");
	*visuclass=0;
        /* Creation d'une nouvelle resource X de type colormap aupres du serveur X */
        /* toutes les cellules de la colormap doivent etre allouees et reservees au 
           client X c-a-d qu'il est le seul a en disposer tant en lecture qu'en ecriture */
        color_map = XCreateColormap ( display_mef, fenetre_mef0, visual, AllocAll );

        /* le nombre de couleurs de la palette colormap */
        *nbcolo = visual->map_entries;  
        if (DefaultDepth(display_mef, screen_mef) == 1) 
        { 
         /*  Black and white: tout ce qui n'est pas noir est blanc */
         /* La palette est blanche */
          for (i = 1; i < *nbcolo; i++)
          { /* au dela 0 toutes les couleurs sont blanches */
            red[i]  =1.0;
            green[i]=1.0;
            blue[i] =1.0;
          }
        }
        else
        { /* Il y a nbcolo>2 couleurs dans la palette
             Remplissage de red green blue a partir de la palette par defaut 
             qui est partageable donc necessairement en lecture seule */ 
          color_map0 = DefaultColormap(display_mef, screen_mef); 
          proc(xvColormapToRGB)( color_map0,  red, green, blue, *nbcolo);
        }  

        /* stockage des red green blue dans color_map la palette courante Read/Write */
        n1core = (int) 0 ; /* Premiere couleur reservee (non modifiee) */
                           /* Derniere couleur reservee = derniere couleur elementaire */
        n1coel = (int) 0 ; /* Premiere couleur elementaire pour le noir */
                           /* Ici aucune couleur reservee a part les couleurs elementaires*/

        /* Definition des couleurs imposees pour toutes les palettes */
        proc(xvCouleursImposees)( n1coel, &ndcoel, red, green, blue );

        /* Stockage de ces couleurs dans la palette colormap */
        proc(xvStockeRGBtoColormap)(*nbcolo, red, green, blue,  norgb, color_map);
        foreground = (unsigned long) ndcoel;     /* blanc */
        background = (unsigned long) n1coel;     /* noir  */
              
        /* le numero des couleurs reservees 0 a ndcore et modifiables ndcore + 1 a nbcolo-1 */
        ndcore = ndcoel ;
        n1coul = ndcore + (int) 1 ;
        ndcoul = (int) *nbcolo-1 ; 
        break;
      } 

      case GrayScale :  /* Classe des ecrans gris avec 8 plans => 256 niveaux de gris  */
      /*--------------     Les palettes colormap sont en Lecture/ecriture non partageables */
      /*                   Le no couleur mefisto = celui de la palette physique            */
      /*                   Tout pixel de l'ecran a un no de couleur dans la palette        */
      /*                   auquel correspond dans la palette 3 valeurs rouge vert bleu tracees */
      /*                   Par exemple si rouge vert bleu sont connus alors */
      /*                   intensite du gris = = 0.3*rouge + 0.59*vert + 0.11*bleu */
      /*                   donne une bonne valeur                                  */
      { 
	if( langage )
	  printf(" GrayScale: NOT TESTED in Mefisto\n");
	else
	  printf(" GrayScale : NON TESTE dans Mefisto\n");

	*visuclass=1;
        color_map  = DefaultColormap( display_mef, screen_mef );
        foreground = WhitePixel ( display_mef, screen_mef );
        background = BlackPixel ( display_mef, screen_mef );
        *nbcolo = visual->map_entries;  /*  *nbcolo    = 256;  */ 

	/* stockage des red green blue dans color_map la palette courante Read/Write */
	n1core = (int) 0 ; /* Premiere couleur reservee (non modifiee) */
        /* Derniere couleur reservee = derniere couleur elementaire */
        n1coel = (int) 0 ; /* Premiere couleur elementaire pour le noir */
        /* Ici aucune couleur reservee a part les couleurs elementaires*/
        /* Definition des couleurs imposees pour toutes les palettes */
        proc(xvCouleursImposees)( n1coel, &ndcoel, red, green, blue );
       
	/* le numero des couleurs reservees et modifiables */
        ndcoel = (int) 1 ;  
        ndcore = (int) 15 ;
        n1coul = (int) 16 ;
        ndcoul = (int) 255 ;
        break;
      } 

      case DirectColor :  /* Classe des ecrans avec >8 plans ecrans couleurs => >256 couleurs  */
      /*---------------      Les palettes colormap sont en Lecture/ecriture non partageables */
      /*                     Tout pixel de l'ecran a un no de couleur dans la palette physique   */
      /*                     qui se decompose en 3 valeurs rouge vert bleu qui */
      /*                     donnent le RVB trace par une fonction du type */
      /*                     no palette= rouge + vert*2**nbplansrouge + bleu*2**(nbplansrouge+vert) */
      {
	*visuclass=2;
	if( langage )
	{
	  printf("DirectColor : NOT TESTED in Mefisto\n");
	  printf("Try to replace DirectColor by PseudoColor...\n");
	}
	else
	{
	  printf("DirectColor : NON teste dans Mefisto\n");
	  printf("Essai de remplacer DirectColor par PseudoColor...\n");
	}
        visual->class=PseudoColor;
        goto relance;
      }

      /*---------------------------------------------------------------------------*/
      /*  Les 3 CLASSES avec une PALETTE NON MODIFIABLE (C'est celle du serveur X) */
      /*---------------------------------------------------------------------------*/

      case  StaticColor: /* idem PseudoColor: Classe ecrans avec 8 plans ecrans couleurs=>256 couleurs */
      /*---------------     La palette colormap est en Lecture seule, partageable, non modifiable */
      /*                    Le no couleur mefisto = celui de la palette physique            */
      /*                    Tout pixel de l'ecran a un no de couleur dans la palette        */
      /*                    auquel correspond dans la palette 3 valeurs rouge vert bleu tracees */
      {
	if( langage )
	  printf("StaticColor : NOT TESTED in Mefisto\n");
	else
	  printf("StaticColor : NON TESTE dans Mefisto\n");
	*visuclass=3;
        color_map = DefaultColormap( display_mef, screen_mef );
        foreground = WhitePixel ( display_mef, screen_mef );
        background = BlackPixel ( display_mef, screen_mef );
        *nbcolo    = visual->map_entries;
        proc(xvColormapToRGB)( color_map,  red, green, blue, *nbcolo);
        /* le numero des couleurs reservees et modifiables */
        n1core = (int) 0 ;
        n1coel = (int) 0 ;
        ndcoel = (int) 1 ;  
        ndcore = (int) 1 ;
        n1coul = (int) 1 ;
        ndcoul = (int) 1 ;
        break;
      }   

      case StaticGray :   /* idem GrayScale mais avec une colormap fixe       */
      /*---------------      La palette colormap est en Lecture seule, partageable, non modifiable */
      /*                     Le no couleur mefisto = celui de la palette physique            */
      /*                     Tout pixel de l'ecran a un no de gris dans la palette => trace  */
      {
        if( langage )
	  printf("StaticGray : NOT TESTED in Mefisto\n");
	else
	  printf("StaticGray : NON teste dans Mefisto\n");

	*visuclass = 4;
        color_map  = DefaultColormap( display_mef, screen_mef );
        foreground = WhitePixel ( display_mef, screen_mef );
        background = BlackPixel ( display_mef, screen_mef );
        *nbcolo    = 2;
        /* le numero des couleurs reservees et modifiables */
        n1core = (int) 0 ;
        n1coel = (int) 0 ;
        ndcoel = (int) 1 ;   
        ndcore = (int) 1 ;
        n1coul = (int) 1 ;
        ndcoul = (int) 1 ;
        break;
      }  

      case TrueColor  : /* Classe des ecrans avec >8 plans ecrans couleurs => >256 couleurs  */
      /*---------------    La palette colormap est en Lecture seule partageable */
      /*                   Le no couleur mefisto est different de celui de la palette physique */
      /*                   Tout pixel de l'ecran a un no de couleur dans la palette physique   */
      /*                   ce no se decompose en 3 valeurs rouge vert bleu qui */
      /*                   donnent le RVB trace par une fonction du type */
      /*                   no palette= rouge + vert*2**nbplansrouge + bleu*2**(nbplansrouge+vert) */
      { 
	if( langage )
	  printf(" Class TrueColor\n");
	else
	  printf(" Classe TrueColor\n");
	*visuclass = 5;
        /* Palette de couleurs en lecture seule donc partageable */
        /* Recuperation de la  colormap du serveur X */
        color_map = DefaultColormap( display_mef, screen_mef );

        /* le nombre de couleurs de la palette mefisto est celui du stockage de red,...*/
        /* Attention: nbcolo n'est pas le nombre de couleurs de la palette physique! */
        *nbcolo = visual->map_entries;
        if (*nbcolo>CMAPSIZE) *nbcolo = (int) CMAPSIZE;
        *nbcolo = 256;
	if( langage )
	  printf(" TrueColor: Mefisto with %d colors\n", *nbcolo);
	else
	  printf(" TrueColor: Mefisto avec %d couleurs\n", *nbcolo);

        /* Il y a nbcolo>2 couleurs dans la palette mefisto
           Remplissage de red green blue a partir de la palette par defaut 
           qui est partageable donc necessairement en lecture seule */  

        /* Initialisation des red green blue du C */
        n1core = (int) 0 ; /* Premiere couleur reservee (non modifiee) */
                           /* Derniere couleur reservee = derniere couleur elementaire */
        n1coel = (int) 0 ; /* Premiere couleur elementaire pour le noir */
                           /* Ici aucune couleur reservee a part les couleurs elementaires*/

        /* Definition des couleurs imposees */
        proc(xvCouleursImposees)( n1coel, &ndcoel, red, green, blue );

        /* le numero des couleurs reservees 0 a ndcore et modifiables ndcore+1 a nbcolo-1 */
        ndcore = ndcoel ;
        n1coul = ndcore + (int) 1 ;
        ndcoul = (int) *nbcolo-1 ;

        /* les couleurs mefisto modifiables sont momentanement grises */
        for (i = n1coul; i <= ndcoul; i++)
        { /* au dela de ndcore toutes les couleurs sont grises */
          red[i]  = (float) i/ndcoul;
          green[i]= (float) i/ndcoul;
          blue[i] = (float) i/ndcoul;
        }

        /* les nbcolo couleurs sont allouees dans la colormap partageable */
        /* le norgb[i] de chaque couleur i de la palette mefisto est celui de la couleur
           la plus approchante dans la palette physique */
        proc(xvStockeRGBtoColormap)(*nbcolo, red, green, blue,  norgb, color_map);

        foreground = (unsigned long) ndcoel;     /* blanc */
        background = (unsigned long) n1coel;     /* noir  */
        break;
      } 
    }

    /* Traitement quelquesoit la classe */
    /*----------------------------------*/
    printf(" Visual->map_entries = %d\n", visual->map_entries);
    if( langage )
    {
      printf(" Used  colors number = %d\n", *nbcolo);
      printf(" Color Planes Depth  = %d\n", DefaultDepth(display_mef,screen_mef) );
    }
    else
    {
      printf(" Nb couleurs utilisees         = %d\n", *nbcolo);
      printf(" Nb de plans pour les couleurs = %d\n", DefaultDepth(display_mef,screen_mef) );
    }


    attributes.background_pixel = (unsigned long) BlackPixel(display_mef,screen_mef);
    attributes.border_pixel     = (unsigned long) WhitePixel(display_mef,screen_mef);
    attributes.backing_store    = Always;
    attributes.colormap         = color_map;

    /*==========================*/
    /* Creation de la fenetre X */
    /*==========================*/
    fenetre_mef = XCreateWindow( display_mef,
			    fenetre_mef0,
                            x, y, width, height, border_width,
                            CopyFromParent,
                            InputOutput,
                            visual,
                            CWBackPixel | CWBorderPixel | CWBackingStore | CWSaveUnder | CWColormap,
                            &attributes);

    strcpy(data,"Mefisto");       /* le nouveau titre Mefisto */
    XChangeProperty( display_mef, fenetre_mef, XA_WM_NAME, XA_STRING, 8,
                     PropModeReplace,(unsigned char *) data, strlen("Mefisto") );

    /* les seuls EVENEMENTS qui interessent MEFISTO  FAITS AU DESSUS */
    /* 18/ 1/1999 attributes.event_mask =  ExposureMask | KeyPressMask | ButtonPressMask ; */
    /* 10/10/2000 ajout de | ButtonMotionMask | ButtonReleaseMask */
    /* 19/01/2001 ajout de | PointerMotionMask  pour tout deplacement du pointeur de la souris */
    /*                                          sans bouton enfonce ou relache */
    attributes.event_mask = VisibilityChangeMask | StructureNotifyMask 
                           | ExposureMask        | KeyPressMask
                           | ButtonPressMask     | ButtonReleaseMask 
                           | ButtonMotionMask    | PointerMotionMask;

    XChangeWindowAttributes( display_mef, fenetre_mef,
                             CWEventMask, &attributes); 

    XMapWindow(display_mef, fenetre_mef); /* la fenetre devient visible */ 

    /* retourne dans les variables FORTRAN les numeros: premiere et derniere couleur
       reservee, numeros premiere et derniere couleur elementaire 
       et couleurs modifiables dans la palette des couleurs mefisto */  
    *n1coref = n1core ;
    *n1coelf = n1coel ;
    *ndcoelf = ndcoel ;
    *ndcoref = ndcore ;  
    *n1coulf = n1coul ;
    *ndcoulf = ndcoul ; 

    /*  Recuperation des fontes X11  */
    /*-------------------------------*/
    /* Chargement des fontes X11 avec le filtre "*"  */
    listfonts = XListFonts (display_mef, "*", *maxfonts, nbfonts);
    /* si necessaire, on prend un filtre plus large */
    if (*nbfonts == 0) listfonts = XListFonts(display_mef, "*-*", *maxfonts, nbfonts);
    if (*nbfonts != 0)
    {
      font_id = XLoadFont (display_mef, *listfonts);
      nbpolices = *nbfonts;
      if( langage )
	printf(" Available X11-Fonts number = %d\n", nbpolices);
      else
	printf(" Nb de FONTES X11 disponibles = %d\n", nbpolices);
    }
    else
    {
      if(langage)
	printf ("NO AVAILABLE X11 FONTS! => Modify X11 installation\n"); 
      else
	printf ("AUCUNE FONTE X11 DISPONIBLE!=>  Revoir installation X11\n"); 
      exit(1);
    }

    gcvalues.foreground = norgb[foreground];
    gcvalues.background = norgb[background];
    gcvalues.function   = GXcopy    ;
    gcvalues.font       = font_id   ;
    gc_mef = XCreateGC( display_mef,
		        fenetre_mef,
		        GCForeground | GCBackground | GCFunction | GCFont, 
		        &gcvalues );
    XSetFillRule( display_mef,gc_mef,WindingRule );

    XSetArcMode( display_mef,gc_mef,ArcChord );

    /*  l'analyse des noms sera ensuite faite en fortran  */ 
    /*  le calcul du nombre de caracteres du nom de chaque fonte */
    for ( k = 0; k < nbpolices; k++ )
    { 
      nbchar[k] = strlen( listfonts[k] );
      strcpy( namefonts[k], listfonts[k] );
  /*  if( langage )
	printf("xvinfo Fonte X11 %d: %s\n", k, namefonts[k]);
      else
        printf("xvinfo X11-Font %d: %s\n", k, namefonts[k]);  */
      }
    
    /* Creation du pixmap ou sont effectues les traces avant Xcopy dans fenetre_mef */
    /*------------------------------------------------------------------------------*/
    mempx = XCreatePixmap( display_mef, fenetre_mef, lapxfe, lhpxfe,
  		           DefaultDepth( display_mef, screen_mef ) );

    /* Creation du pixmap pour sauvegarde restauration de la fenetre_mef de Mefisto */
    /*------------------------------------------------------------------------------*/
    mempxsauvfen = XCreatePixmap( display_mef, fenetre_mef, lapxfe, lhpxfe,
  		        	  DefaultDepth( display_mef, screen_mef ) );
    
    /* Creation du mempxaccro pour le carre d'accrochage d'un item  */
    /*--------------------------------------------------------------*/
    proc(initaccrochage)();
}

void proc(xvrecuprgbdec)(int *nbcolor, float *r, float *g, float *b)  
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :  INITIALISE les tableaux FORTRAN r,v,b, a partir des tableaux 
 -----  red, green, blue, de la palette XVUE
      ( Transfert de la palette C des couleurs X dans la palette FORTRAN )

 ENTREE :
 --------
 nbcolor : nombre de couleurs dans la palette des couleurs a recuperer en FORTRAN

 SORTIES :
 ---------
 r, g, b : les tableaux ROUGE VERT BLEU de la palette en FORTRAN
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/
{
  int i;

  for (i=0; i< *nbcolor; i++)
  {
    *(r+i) = red[i]  ;
    *(g+i) = green[i];
    *(b+i) = blue[i] ;
  }
}


void proc(xvactivervb)( int *palcour, int *nbcells,
                        float r[], float g[], float b[] )   
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :     ACTIVER dans la colormap physique les nbcells couleurs definies 
 -----     dans r, g, b de FORTRAN  et les recopier dans red green blue de C
           ( cf son appel dans ~/xvue/palcde.f )

 ENTREES :
 ---------
 palcour : numero de la palette courante
 nbcells : nombre de couleurs dans la palette des couleurs a charger
 r, g, b : les tableaux ROUGE VERT BLEU (parametres FORTRAN) de la palette

 SORTIES dans les TABLEAUX C GLOBAUX de xvue:
 --------------------------------------------
 norgb   : numero des couleurs de la palette mefisto dans la palette physique
           tableau interne au C de mefisto
 red, green, blue : les tableaux ROUGE VERT BLEU des tableaux C
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                 AOUT 1997
12345X7..............................................................012345678*/ 
{
  int i;

  /* le numero de la palette courante est palcour */
  palcourc = *palcour ;

  /* les r,g,b de >FORTRAN sont installes dans la colormap 
     en retour, initialisation du tableau C norgb
     norgb[i]=numero de la couleur i dans la palette mefisto 
              dans la palette physique  */
  proc(xvStockeRGBtoColormap)( *nbcells, r, g, b, 
                               norgb, color_map );

  XInstallColormap( display_mef, color_map );
  XFlush( display_mef );

  /* Transfert de r,v,b de FORTRAN dans les tableaux red green blue de C */
  for ( i = 0 ; i < *nbcells ; i++ )
  {
    red[i]   = r[i];
    green[i] = g[i];
    blue[i]  = b[i];
  }
}


void proc(xvcouleur)(int *icolor)
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :     IMPOSER pour tous les traces a venir la couleur icolor de la palette
 -----     ( appelable en FORTRAN )
            
 ENTREES :
 ---------
 icolor  :  le numero dans la palette mefisto de la couleur de trace
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                 AOUT 1997
12345X7..............................................................012345678*/ 
{  
  unsigned long  colour;

  /* indirection necessaire pour traiter a la fois PseudoColor et TrueColor */
  /* colour est le numero de la couleur icolor mefisto dans la palette physique */
  if (*icolor < 0) { *icolor=1; } /* rouge */
  if (*icolor >= CMAPSIZE ) { *icolor=1; } /* rouge */
  colour = norgb[*icolor];
  foreground = *icolor;
  XSetForeground (display_mef, gc_mef, colour);

  /* Traitement du postscript :
     counb est ramenee sur [0,1] suivant la palette utilisee
     counb = -1 si *icolor est une des couleurs non modifiable de la palette */
  if (lasopsc > 0)
    {
      if (icolorm != *icolor)
	{ 
	  icolorm = *icolor ;  
	  if (nbrcon > 0)
	    {   
	      fprintf(fpo,"%s",concat) ;
	      nbrcon  = 0 ;
	      concat[0] = '\0' ;
	    } 
	  courgb[0] = red[*icolor] ;
	  courgb[1] = green[*icolor] ;
	  courgb[2] = blue[*icolor] ;
	  if (*icolor < n1coul)
	    {
	      counb = -1 ;
	    }
	  else
	    {
	      if (palcourc == 2)
		{
		  /* *icolor varie sur palette de 6 couleurs */
		  counb = ( (*icolor-n1coul) % 6 ) / 5.0  ;
		}
	      else
		{
		  if (palcourc == 12)
		    {
		      /* *icolor varie sur palette de 10 couleurs */
		      counb = ( (*icolor-n1coul) % 10 ) / 9.0 ;
		    }
		  else
		    {
		      /* *icolor varie sur palette totale */
		      counb = ((float) (ndcoul - *icolor)) / (ndcoul - n1coul) ;
		    }
		}
	    }
	}
    }
}

void proc(xvpostscript)( int *lasops )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :   IMPOSER la nouvelle option de sortie du trace du maillage en POSTSCRIPT
 -----       ( appelable en FORTRAN )
         A NE PAS APPELER DIRECTEMENT LORSQU'ON UTILISE xvue
           ==========================
 ENTREE  :
 ---------  
 lapops  : DEMANDE OU NON DE SORTIE DU TRACE DU MAILLAGE EN POSTSCRIPT
           0     : PAS DE SORTIE POSTSCRIPT
           NON 0 : SORTIE DU TRACE EN POSTSCRIPT   
           1 trace normal
           2 trace normal avec menu
           3 ..  10 trace dans les differents menus
          -3 .. -10 effacement des differents menus
           11  12   trace menu qualites ou histos (resp 1 ou 2)
          -11 -12   effacement menu qualites ou histos (resp 1 ou 2)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : DOURSAT CHRISTOPHE UPMC ANALYSE NUMERIQUE PARIS              JUIN 1994
12345X7..............................................................012345678*/ 
{       
  int i ;
  /*int j;*/         
  /*char nom[3] ;*/

  /* Traitement du postscript */

  /* ouverture du postscript */
  if ( *lasops > 0 && lasopsc == 0 ) {
    lasopsc = *lasops ;
    if (fpo != NULL) {        
      fclose(fpo) ;
      fpo = NULL ; 
    }
    if ( modepsc != 0 ){
      if ((fpo = fopen("TEMPORAIRE.QUA","r"))!=NULL) {
        fclose(fpo) ;
        fpo = NULL ; 
        remove("TEMPORAIRE.QUA") ;
      }
    }
    if ((fpo = fopen("TEMPORAIRE.EPS","w"))==NULL)
    {
      if( langage )
	printf ("Error: NO creation of file... IMMEDIATE EXIT\n") ;
      else
	printf ("Erreur: NON creation de fichier... SORTIE IMMEDIATE\n") ;
      exit (1);
    }
    if ( modepsc != 0 ){
      for (i = 0; i < 8 ; i++) chaine[i] = (char*) calloc(longchaine[i],sizeof(char)) ;
    }
    nbrcon  = 0 ;
    concat[0] = '\0' ;
    iTe = 0 ;  iFa = 0 ;  ity = 0 ;  iep = 0 ;  iPo = 0 ;
    ire = 0 ;  iRe = 0 ;  iel = 0 ;  iEl = 0 ;  iFP = 0 ;
  }
  else {
    if ( *lasops == 0 ) {
      lasopsc = *lasops ;
      if (fpo != NULL) {        
        fprintf(fpo,"%s",concat) ;
        fclose(fpo) ;
        fpo = NULL ; 
        if ( modepsc != 0 ){
          for (i = 0; i < 8 ; i++) { chaine[i] = '\0' ; free(chaine[i]) ; }
        }
        menu = 0 ;
      }
    }
    else {
      if (*lasops > 100) {
        /* remise a zero des controles de macros */
        iTe = 0 ;  iFa = 0 ;  ity = 0 ;  iep = 0 ;  iPo = 0 ;
        ire = 0 ;  iRe = 0 ;  iel = 0 ;  iEl = 0 ;  iFP = 0 ;
        /* effacement des differents fichiers temporaires */
        if (fpo != NULL) {        
          fclose(fpo) ;
          fpo = NULL ; 
        }
        /* effacement du fichier temporaire des histogrammes de qualite */
        if ( modepsc != 0 ){
          if ((fpo = fopen("TEMPORAIRE.QUA","r"))!=NULL) {
            fclose(fpo) ;
            fpo = NULL ; 
            remove("TEMPORAIRE.QUA") ;
          }
        }
        /* effacement du fichier temporaire postscript */
        if ((fpo = fopen("TEMPORAIRE.EPS","w"))==NULL)
	{
	  if( langage )
	    printf ("Error NO creation of file... EXIT immediately\n") ;
	  else
	    printf ("Erreur en creation de fichier...Sortie immediate\n") ;
	  exit (1);
        }
        /* effacement des menus postscript */
        if ( modepsc != 0 ){for (i = 0; i < 8 ; i++) *chaine[i] = '\0' ;}
        lasopsc = lasopsc - 100 ;
      }
      else {
        fprintf(fpo,"%s",concat) ;
        nbrcon  = 0 ;
        concat[0] = '\0' ;
        if (*lasops < -1) {
          /* effacement du menu correspondant */
          lasopsc = max(*lasops,-11) ;
          if ( modepsc != 0 ){*chaine[-lasopsc-4] = '\0' ;}
        }
        else{
          lasopsc = min(*lasops,11) ;
          if (lasopsc == 2) menu = 1 ;
        }
      }
    }
  }
} 


void proc(fenetremempx)()
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :     Copier la pixmap fenetre_mef dans la pixmap mempx
 -----     ( appelable en FORTRAN )
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN Labo ANALYSE NUMERIQUE LJLL UPMC PARIS   NOVEMBRE 2003
12345X7..............................................................012345678*/ 
{			  
  XSetFunction( display_mef, gc_mef, GXcopy );
  XCopyArea( display_mef, fenetre_mef, mempx, gc_mef,
             0, 0, lapxfe, lhpxfe, /* source */
	     0, 0 );               /* destination */
}

void proc(mempxfenetre)()
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :     Copier la pixmap mempx dans la pixmap fenetre_mef
 -----     ( appelable en FORTRAN )
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN Labo ANALYSE NUMERIQUE LJLL UPMC PARIS   NOVEMBRE 2003
12345X7..............................................................012345678*/ 
{		  
  XSetFunction( display_mef, gc_mef, GXcopy );
  XCopyArea( display_mef, mempx, fenetre_mef, gc_mef,
	     0, 0, lapxfe, lhpxfe, /* source */
	     0, 0 );               /* destination */
}

void proc(sauvefenetre)()
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :     Sauvegarder dans un pixmap les pixels de la fenetre Mefisto
 -----     ( appelable en FORTRAN )
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS              OCTOBRE 2000
12345X7..............................................................012345678*/ 
{
  /* Sauvegarde de la fenetre Mefisto dans mempxsauvfen */			  
  XSetFunction( display_mef, gc_mef, GXcopy );
  XCopyArea( display_mef, fenetre_mef, mempxsauvfen, gc_mef,
             0, 0, lapxfe, lhpxfe, /* source */
	     0, 0 );               /* destination */
}

void proc(restaurefenetre)()
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :     Restaurer les pixels de la fenetre Mefisto a partir d'un pixmap
 -----     ( appelable en FORTRAN )
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS              OCTOBRE 2000
12345X7..............................................................012345678*/ 
{
  /* Restauration de mempxsauvfen dans la fenetre Mefisto */			  
  XSetFunction( display_mef, gc_mef, GXcopy );
  XCopyArea( display_mef, mempxsauvfen, fenetre_mef, gc_mef,
	     0, 0, lapxfe, lhpxfe, /* source */
	     0, 0 );               /* destination */
}

void proc(sauvemempx)()
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :     Sauvegarder dans un pixmap les pixels de mempx
 -----     ( appelable en FORTRAN )
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS              OCTOBRE 2000
12345X7..............................................................012345678*/ 
{
  /* Sauvegarde de la fenetre Mefisto dans mempxsauvfen */			  
  XSetFunction( display_mef, gc_mef, GXcopy );
  XCopyArea( display_mef, mempx, mempxsauvfen, gc_mef,
             0, 0, lapxfe, lhpxfe, /* source */
	     0, 0 );               /* destination */
}

void proc(restauremempx)()
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :     Restaurer les pixels de mempx
 -----     ( appelable en FORTRAN )
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS              OCTOBRE 2000
12345X7..............................................................012345678*/ 
{
  /* Restauration de mempxsauvfen dans la fenetre Mefisto */			  
  XSetFunction( display_mef, gc_mef, GXcopy );
  XCopyArea( display_mef, mempxsauvfen, mempx, gc_mef,
	     0, 0, lapxfe, lhpxfe, /* source */
	     0, 0 );               /* destination */
}

void proc(effacemempx)()
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
but: Effacer le contenu de mempx (mais pas la fenetre_mef)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
{
  /*effacement de la MemoirePx avec la couleur de fond*/
  XSetForeground( display_mef, gc_mef, norgb[background] );
  XFillRectangle( display_mef, mempx, gc_mef, 0, 0, lapxfe, lhpxfe );

  /* Traitement du postscript :
     101 ou 102 effacement en mode sans (1) ou avec (2) menus */
  if (lasopsc > 0)
  {
    lasopsc = 100 + lasopsc ;
    proc(xvpostscript)(&lasopsc) ;
  }
}

void proc(effacer)() 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
 BUT :     EFFACER la fenetre actuelle avec la couleur de fond
 -----  
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  XClearWindow( display_mef, fenetre_mef );
  XFlush( display_mef );
  proc(fenetremempx)();

  /* Traitement du postscript :
     101 ou 102 effacement en mode sans (1) ou avec (2) menus */
  if (lasopsc > 0)
  {
    lasopsc = 100 + lasopsc ;
    proc(xvpostscript)(&lasopsc) ;
  }
}
 
void proc(xvfond)(int *icolor)  
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :     IMPOSER une couleur comme FOND de FENETRE
 -----     ( appelable en FORTRAN )
            
 ENTREES :
 ---------
 icolor  :  le numero dans la palette mefisto de la couleur du fond
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                 AOUT 1997
12345X7..............................................................012345678*/ 
{  
  unsigned long  colour;

  /* indirection necessaire pour traiter a la fois PseudoColor et TrueColor */
  /* colour est le numero de la couleur icolor mefisto dans la palette physique */ 
  background =       *icolor;
  colour     = norgb[*icolor];

  attributes.background_pixel = colour;
  XChangeWindowAttributes( display_mef, fenetre_mef,
                           CWBackPixel | CWBorderPixel | CWBackingStore | CWColormap,
                           &attributes); 

  XSetBackground ( display_mef, gc_mef, colour );
  XSetWindowBackground ( display_mef, fenetre_mef, colour );
}


void proc(xvchargefonte)( int *nofont0, int *nofont, int *largpx, int *hautpx )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :     CHARGER la fonte n de nom NameFont  parmi les fontes disponibles
 -----     ( Appelable en FORTRAN )

 ENTREES :
 ---------
 nofont0 : numero de la fonte actuelle a decharger
 nofont  : numero de la fonte a charger

 SORTIES :
 ---------   
 largpx  : nombre de pixels en largeur de la fonte
 hautpx  : nombre de pixels en hauteur de la fonte
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
 MODIF  : PERRONNET ALAIN LJLL UPMC PARIS & Saint PIERRE DU PERRAY  OCTOBRE 2013
12345X7..............................................................012345678*/ 
{ 
  int    i, ideb, ifin, ilon ;
  char  *font, mp[4] ;    /* Modif:  mp[4] au lieu de mp[3] 10/1/2010 */
  int    direction ;             

  if ( ( 0 <= *nofont ) && (*nofont <= nbpolices - 1 ) )
  {
    /* Traitement de la fonte X11 */
    if ( *nofont0 > 0 )
    {
      /* printf("xvchargefonte decharge la fonte %d\n", *nofont0); */
      XFreeFont( display_mef, struc_police );
      /*  XUnloadFont( display_mef, struc_police->fid ); */
    }

    /* printf("xvchargefonte   charge la fonte %d\n", *nofont); */
    font = listfonts[*nofont];
    struc_police = XLoadQueryFont( display_mef, font );
    if (struc_police != NULL) 
      XSetFont( display_mef, gc_mef, struc_police->fid ) ;

    sprintf(&format[0],"abcdefghijklmnopqrstuvwxyz ABCDEFGHIJKLMNOPQRSTUVWXYZ 1234567890") ;
    XTextExtents( struc_police , format , strlen(format) , &direction ,
                  &ascent_pol , &descent_pol , &mesure ) ;
    *largpx = mesure.width / 64.0 ;
    /* *largpx = *largpx + 1 ; */
    *hautpx = ascent_pol + descent_pol ;  /* VALEUR MAX TROP GRANDE  NON EXACTE  Mai 2012*/
    /* *hautpx = mesure.ascent + mesure.descent ; */

    /* Traitement de la fonte equivalente en PostScript */
    if (lasopsc > 0)
    {
      ideb = 0 ; ifin = 0 ;
      for ( i = 1 ; i < strlen(font) ; i++ )
	{ if (font[i] == '-')
	    { if (ifin == 0) { ideb = i ; ifin = i ; }
	      else
		{ ifin = i ; sprintf(format,"/%%.%ds ",ifin-ideb-1) ;       
		  sprintf(fontcour,format,&font[ideb+1]) ;
		  if (strcmp("/new century schoolbook ",fontcour) == 0 )
		    { sprintf(fontcour,"/NewCenturySchlbk ") ; }
		  fontcour[1] = toupper( fontcour[1]) ; break ; } } }
      ideb = ifin ;
      for ( i = ideb + 1 ; i < strlen(font) ; i++ )
	{ if (font[i] == '-') { ifin = i ; sprintf(format,"/%%.%ds ",ifin-ideb-1) ;       
	    ilon = strlen(fontcour) ; sprintf(&fontcour[ilon],format,&font[ideb+1]) ;       
	    fontcour[ilon+1] = toupper(fontcour[ilon+1]) ; break ; } }
      ideb = ifin ; 
      for ( i = ideb + 1 ; i < strlen(font) ; i++ )
	{ if (font[i] == '-') { ifin = i ; sprintf(format,"%.1s",&font[ideb+1]) ;       
	    if (format[0] == 'o') { sprintf(&fontcour[strlen(fontcour)],"/Oblique ") ;
	    } else { if (format[0] == 'i') { sprintf(&fontcour[strlen(fontcour)],"/Italic ") ;
	      } else { if (format[0] == 'r') { sprintf(&fontcour[strlen(fontcour)],"/ ") ;
		       } else { sprintf(&fontcour[strlen(fontcour)],"/  ") ; } } } break ; } }
      ideb = ifin ; ilon = 0 ;
      for ( i = ideb + 1 ; i < strlen(font) ; i++ )
	{ if (font[i] == '-')
	    { ilon = ilon + 1 ;
	    if (ilon < 7) { ideb = i ; ifin = i ;
	    }
	    else
	    { ifin = i ;
	      sprintf(format,"(%%.%ds)",ifin-ideb-1) ;     
	      sprintf(mp,format,&font[ideb+1]) ; 
	      if (strcmp("(p)",mp) != 0 ) {
		sprintf(mp,"(m)") ; }
	      break ; } } }
      sprintf(&fontcour[strlen(fontcour)],"%d %d %s charge\n",mesure.ascent,
	      mesure.rbearing-mesure.lbearing,mp) ; 
      
      if (nbrcon > 0)
      {   
	fprintf(fpo,"%s",concat) ;
	nbrcon  = 0 ;
	concat[0] = '\0' ;
      }
      if (lasopsc < 3)
      {
	fprintf(fpo,"%s",fontcour) ;
      }
      else
      {
	sprintf(&chaine[lasopsc-4][strlen(chaine[lasopsc-4])],"%s",fontcour) ;
      }
    }
  }    
  else
  {
    if( langage )
      printf("\n xvchargefonte: FONTE INDISPONIBLE %d => FONTE INCHANGEE\n",*nofont);
    else
      printf("\n xvchargefonte: %d UNAVAILABLE FONT => FONT UNCHANGED\n",*nofont);
  }
} 

void proc(xvnbpixeltexte)( char *texte, int *length, int *nbpxla, int *nbpxha )  
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :     Retrouver les dimensions en pixels d'une chaine de caracteres
 -----     ( Appelable en FORTRAN )
            
 ENTREES :
 ---------  
 texte   : la chaine de caracteres 
 length  : nombre de caracteres de texte 

 SORTIES :
 ---------
 nbpxla  : le nombre de pixels de la largeur de la chaine de caracteres
 nbpxha  : le nombre de pixels de la hauteur de la chaine de caracteres
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : DOURSAT CHRISTOPHE UPMC ANALYSE NUMERIQUE PARIS          NOVEMBRE 1994
12345X7..............................................................012345678*/ 
{        
  int  direction ; 

  XTextExtents(struc_police , texte , *length , &direction ,
               &ascent_pol , &descent_pol , &mesure ) ;  
  *nbpxla = mesure.width ;
  *nbpxha = mesure.ascent + mesure.descent ;
}

void proc(xvfermer)()   
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :     FERMER toutes les ressources X11 utilisees par Mefisto
 -----  
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
 MODIFS : PERRONNET ALAIN LJLL UPMC PARIS & Saintt PIERRE DU PERRAY OCTOBRE 2013
12345X7..............................................................012345678*/ 
{
 /* XUnmapWindow( display_mef, fenetre_mef); la fenetre devient invisible */ 
  XFreeFontNames( listfonts );
  XFreeGC(        display_mef, gc_mef );
  XFreeColormap(  display_mef, color_map );
  XDestroyWindow( display_mef, fenetre_mef );
  XCloseDisplay(  display_mef );
}

void proc(xvpxfenetre)( int *x, int *y )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :     RECUPERER le nombre de pixels de la largeur et hauteur de 
 -----     la FENETRE ACTUELLE X  
           ( appelable en FORTRAN )

 SORTIES :
 ---------
 x       : nombre de pixels de la largeur de la fenetre
 y       : nombre de pixels de la hauteur de la fenetre
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  XWindowAttributes wattrib;

  XGetWindowAttributes(display_mef, fenetre_mef, &wattrib); 
  *x = wattrib.width;
  *y = wattrib.height;
  xpixels = *x ;
  ypixels = *y ;
}

void proc(xvftexte)( char string[], int *length, int *x1, int *y1 )   
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :     TRACER un texte a partir du pixel (x1,y1) dans la fenetre_mef
 -----         ( appelable en FORTRAN )

 ENTREES :
 ---------  
 string  : la chaine de caracteres a tracer
 length  : le nombre de caracteres de la chaine
 x1      : le numero du pixel en largeur du debut de la chaine
 y1      : le numero du pixel en hauteur du debut de la chaine 
           (origine au coin superieur gauche de la fenetre)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{   
  XDrawString(display_mef, fenetre_mef, gc_mef, *x1, *y1, string, *length);
 
}

void proc(xvtexte)( char string[], int *length, int *x1, int *y1 )   
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :     TRACER un texte a partir du pixel (x1,y1) de la fenetre
 -----         ( appelable en FORTRAN )

 ENTREES :
 ---------  
 string  : la chaine de caracteres a tracer
 length  : le nombre de caracteres de la chaine
 x1      : le numero du pixel en largeur du debut de la chaine
 y1      : le numero du pixel en hauteur du debut de la chaine 
           (origine au coin superieur gauche de la fenetre)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{   
  XDrawString(display_mef, mempx, gc_mef, *x1, *y1, string, *length);
 /*  proc(xvvoir)(); */

  /* Traitement du postscript */
  if (lasopsc > 0){
    iTe = 1 ;
    if (nbrcon > 0){   
      fprintf(fpo,"%s",concat) ;
      nbrcon  = 0 ;
      concat[0] = '\0' ;
    }
    sprintf(format,"(%%.%ds) %%6i %%6i %%4.2f %%4.2f %%4.2f 0.00 T\n",*length) ;
    if (lasopsc < 3){
      fprintf(fpo,format,string,*x1,ypixels-*y1,courgb[0],courgb[1],courgb[2]) ;
    }
    else{
      sprintf(&chaine[lasopsc-4][strlen(chaine[lasopsc-4])],
              format,string,*x1,ypixels-*y1,courgb[0],courgb[1],courgb[2]) ;
    }                             
  }
}  


void proc(xvface)( int *n, XPoint *pts )  /*  ATTENTION XPoint p => short  p.x et p.y */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
 BUT :     TRACER le remplissage du polygone defini par n points pt
 -----       ( appelable en FORTRAN )

 ENTREES :
 ---------  
 n       : nombre de sommets du polygone
 pts     : coordonnees pixels (INTEGER*2 !) des n sommets du polygone
           le pt 1 est raccorde au pt n
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  int npts, i , j;

  XFillPolygon(display_mef, mempx, gc_mef, pts, *n, Complex, CoordModeOrigin);
  /* proc(xvvoir)(); */

  /* Traitement du postscript :
     instruction ps : F
     destination : TEMPORAIRE.EPS si l'instruction est de type dessin
                 : chaine   si l'instruction correspond a un menu */
  if (lasopsc > 0){ 
    npts = *n;
    iFa = 1 ;
    if (nbrcon > 0){   
      fprintf(fpo,"%s",concat) ;
      nbrcon  = 0 ;
      concat[0] = '\0' ;
    }
    buf[0] = '\0' ;
    for (j = 0 ; j <= npts / 16 ; j++)
    {
      for (i = j * 16 ; i <= min( npts-1 , (j+1) * 16 -1 ) ; i++) 
        sprintf(&buf[strlen(buf)], "%6i %6i " , pts[i].x , ypixels-pts[i].y) ;
      if (j == npts / 16) {
        if (counb != -1) {
          sprintf(&buf[strlen(buf)],"%3i %4.2f %4.2f %4.2f %4.2f F\n",npts,courgb[0],courgb[1],courgb[2],counb) ;
        }
        else {
          sprintf(&buf[strlen(buf)],"%3i %4.2f %4.2f %4.2f 1.00 F\n",npts,courgb[0],courgb[1],courgb[2]) ;
        }
      }
      else {
        sprintf(&buf[strlen(buf)],"\n") ;
      }
      if (lasopsc < 3) {
        fprintf(fpo,"%s",buf) ;
      }
      else {
        sprintf(&chaine[lasopsc-4][strlen(chaine[lasopsc-4])],"%s",buf) ;
      }
      buf[0] = '\0' ;
    }
  }
} 


void proc(xvtypetrait)( int *ptype )   
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
 BUT :     IMPOSER le type de trace d'un trait
 -----        ( appelable en FORTRAN )

 ENTREE :
 --------
 ptype   : 0 LIGNE CONTINUE
           1 LIGNE TIRETEE
           2 LIGNE TIRETEE D'EPAISSEUR DOUBLE
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  XGCValues gcvalues;
  int type;

  type=  *ptype;
  switch (type)
  {
    case 0  : {gcvalues.line_style = LineSolid;     break;}
    case 1  : {gcvalues.line_style = LineOnOffDash; break;}
    default : {gcvalues.line_style = LineDoubleDash;break;}
  }
  XChangeGC(display_mef, gc_mef, GCLineStyle, &gcvalues);

  /* Traitement du postscript :
     destination : TEMPORAIRE.EPS si l'instruction est de type dessin
                 : chaine   si l'instruction correspond a un menu */
  if ( lasopsc > 0){ 
    ity = 1 ;
    if (nbrcon > 0){   
      fprintf(fpo,"%s",concat) ;
      nbrcon  = 0 ;
      concat[0] = '\0' ;
    }
    sprintf(buf,"%2i typet\n",type) ;
    if (lasopsc < 3){
      fprintf(fpo,"%s",buf) ;
    }
    else{
      sprintf(&chaine[lasopsc-4][strlen(chaine[lasopsc-4])],"%s",buf) ;
    }
  }
}  
    

void proc(xvepaisseur)( int *pepais )   
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
 BUT :     IMPOSER le nombre d'epaisseurs des traits
 -----         ( appelable en FORTRAN )
              
 ENTREE :
 --------
 pepais : nombre d'epaisseurs des traces a venir
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  XGCValues gcvalues;

  gcvalues.line_width = *pepais;
  XChangeGC(display_mef, gc_mef, GCLineWidth, &gcvalues);

  /* Traitement du postscript :
     destination : TEMPORAIRE.EPS si l'instruction est de type dessin
                 : chaine   si l'instruction correspond a un menu */
  if ( lasopsc > 0){ 
    iep = 1 ;
    if (nbrcon > 0){   
      fprintf(fpo,"%s",concat) ;
      nbrcon  = 0 ;
      concat[0] = '\0' ;
    }
    if (lasopsc < 3){
      fprintf(fpo,"%2i epais\n" , *pepais) ;
    }
    else{
      sprintf(&chaine[lasopsc-4][strlen(chaine[lasopsc-4])],
              "%2i epais\n" , *pepais) ;
    }
  }
}
  

void proc(xvftrait)( int *x1, int *y1, int *x2, int *y2 ) 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :   TRACER du point (x1,y1) PIXELS au point (x2,y2) PIXELS dans fenetre_mef
 -----         ( appelable en FORTRAN )

 ENTREES :
 ---------
 x1,y1   : coordonnees pixels du point INITIAL (ORIGINE coin superieur gauche)
 x2,y2   : coordonnees pixels du point FINAL   (ORIGINE coin superieur gauche)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{   
  XDrawLine( display_mef, fenetre_mef, gc_mef, *x1, *y1, *x2, *y2 );
}


void proc(xvtrait)( int *x1, int *y1, int *x2, int *y2 ) 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :     TRACER du point (x1,y1) PIXELS au point (x2,y2) PIXELS
 -----         ( appelable en FORTRAN )

 ENTREES :
 ---------
 x1,y1   : coordonnees pixels du point INITIAL (ORIGINE coin superieur gauche)
 x2,y2   : coordonnees pixels du point FINAL   (ORIGINE coin superieur gauche)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{   
  XDrawLine( display_mef, mempx, gc_mef, *x1, *y1, *x2, *y2 ); 
 /*  proc(xvvoir); */

  /* Traitement du postscript :
     instruction ps : S
     destination : TEMPORAIRE.EPS si l'instruction est de type dessin
                 : chaine   si l'instruction correspond a un menu */
  if (lasopsc > 0){ 
    buf[0] = '\0' ;
    if (lasopsc < 3){
      /* On essaye de rendre plus compacts les traces ps de segments consecutifs
         qui sont stockes dans concat. */
      if (nbrcon == 0){     
        nbrcon = 1 ;
        xinic  = *x1 ;
        yinic  = *y1 ;
        xcouc  = *x2 ;
        ycouc  = *y2 ;
        if (counb != -1) {
          sprintf(&concat[0],"%6i %6i %6i %6i %3i %4.2f %4.2f %4.2f %4.2f S\n",
                  *x1, ypixels-*y1, *x2, ypixels-*y2,nbrcon,courgb[0],courgb[1],courgb[2],counb) ;
        }
        else {
          sprintf(&concat[0],"%6i %6i %6i %6i %3i %4.2f %4.2f %4.2f 0.00 S\n",
                  *x1, ypixels-*y1, *x2, ypixels-*y2,nbrcon,courgb[0],courgb[1],courgb[2]) ;
        }
      }
      else {
        if (*x1 == xcouc && *y1 == ycouc){
          nbrcon++ ;
          if (*x2 == xinic && *y2 == yinic){
            iPo = 1 ;
            /* la suite de segments est en fait un contour ferme */
            if (counb != -1) {
              sprintf(&concat[strlen(concat)-26],"%3i %4.2f %4.2f %4.2f %4.2f P\n",
                      nbrcon,courgb[0],courgb[1],courgb[2],counb) ;
            }
            else {
              sprintf(&concat[strlen(concat)-26],"%3i %4.2f %4.2f %4.2f 0.00 P\n",
                      nbrcon,courgb[0],courgb[1],courgb[2]) ;
            }
            /* ecriture de concat */
            fprintf(fpo,"%s",concat) ;
            nbrcon  = 0 ;
            concat[0] = '\0' ;
          } 
          else {
            if (nbrcon % 16 == 0)
            /* concat est deja rempli, on le vide */
            {
              sprintf(&concat[strlen(concat)-26],"\n ") ;
              fprintf(fpo,"%s",concat) ;
              concat[0] = '\0' ;
              sprintf(&concat[0],"                            ") ;
            }
            if (counb != -1) {
              sprintf(&concat[strlen(concat)-26],"%6i %6i %3i %4.2f %4.2f %4.2f %4.2f S\n",
                      *x2, ypixels-*y2,nbrcon,courgb[0],courgb[1],courgb[2],counb) ;
            }
            else {
              sprintf(&concat[strlen(concat)-26],"%6i %6i %3i %4.2f %4.2f %4.2f 0.00 S\n",
                      *x2, ypixels-*y2,nbrcon,courgb[0],courgb[1],courgb[2]) ;
            }
            xcouc  = *x2 ;
            ycouc  = *y2 ;
          }
        }
        else {
          /* ni suite, ni fermeture d'un polygone. c'est donc un nouveau segment */
   
       fprintf(fpo,"%s",concat) ;
          nbrcon = 1 ;
          xinic  = *x1 ;
          yinic  = *y1 ;
          xcouc  = *x2 ;
          ycouc  = *y2 ;
          if (counb != -1) {
            sprintf(&concat[0],"%6i %6i %6i %6i %3i %4.2f %4.2f %4.2f %4.2f S\n",
                   *x1, ypixels-*y1, *x2, ypixels-*y2,nbrcon,courgb[0],courgb[1],courgb[2],counb) ;
          }
          else {
            sprintf(&concat[0],"%6i %6i %6i %6i %3i %4.2f %4.2f %4.2f 0.00 S\n",
                   *x1, ypixels-*y1, *x2, ypixels-*y2,nbrcon,courgb[0],courgb[1],courgb[2]) ;
          }
        }
      } 
    }
    else{
      /* Pour les menus les segments sont traces en ps tels que, un par un */
      if (counb != -1) {
        sprintf(&buf[0],"%6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f S\n",
                *x1, ypixels-*y1, *x2, ypixels-*y2,courgb[0],courgb[1],courgb[2],counb) ;
      }
      else {
        sprintf(&buf[0],"%6i %6i %6i %6i %4.2f %4.2f %4.2f 0.00 S\n",
                *x1, ypixels-*y1, *x2, ypixels-*y2,courgb[0],courgb[1],courgb[2]) ;
      }
      sprintf(&chaine[lasopsc-4][strlen(chaine[lasopsc-4])],"%s",buf) ;
    }
  }
} 

void proc(xvtraits)( int *nbpoints, XPoint *points ) 
/*  ATTENTION XPoint p => short  p.x et p.y */ 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :     TRACER les segments points(i) points(i+1) pour i=1 a nbpoints-1
 -----         ( appelable en FORTRAN )

 ENTREES :
 ---------
 nbpoints : nombre de points sommets des traits
 points   : coordonnees pixels (INTEGER*2 !) des nbpoints points
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{   
  int npts, i , j;

  XDrawLines( display_mef, mempx, gc_mef, points, *nbpoints, CoordModeOrigin );
 /*  proc(xvvoir)();  */

  /* Traitement du postscript :
     instruction ps : P
     destination : TEMPORAIRE.EPS si l'instruction est de type dessin
                 : chaine   si l'instruction correspond a un menu */
  if (lasopsc > 0){ 
    iPo = 1 ;
    npts = *nbpoints-1;
    if (nbrcon > 0){   
      fprintf(fpo,"%s",concat) ;
      nbrcon  = 0 ;
      concat[0] = '\0' ;
    }
    buf[0] = '\0' ;
    for (j = 0 ; j <= npts / 16 ; j++)
    {
      for (i = j * 16 ; i <= min( npts-1 , (j+1) * 16 -1 ) ; i++) 
        sprintf(&buf[strlen(buf)], "%6i %6i " , points[i].x , ypixels-points[i].y) ;
      if (j == npts / 16) {
        if (counb != -1) {
          sprintf(&buf[strlen(buf)],"%3i %4.2f %4.2f %4.2f %4.2f P\n",npts,courgb[0],courgb[1],courgb[2],counb) ;
        }
        else {
          sprintf(&buf[strlen(buf)],"%3i %4.2f %4.2f %4.2f 0.00 P\n",npts,courgb[0],courgb[1],courgb[2]) ;
        }
      }
      else {
        sprintf(&buf[strlen(buf)],"\n") ;
      }
      if (lasopsc < 3) {
        fprintf(fpo,"%s",buf) ;
      }
      else {
        sprintf(&chaine[lasopsc-4][strlen(chaine[lasopsc-4])],"%s",buf) ;
      }
      buf[0] = '\0' ;
    }
  }
}

void proc(xvfacetraits)( int *ncf, int *nca, int *n, XPoint *pts ) 
/*  ATTENTION XPoint p => short  p.x et p.y */ 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
 BUT :     TRACER le remplissage du polygone defini par n points pt
 -----       ( appelable en FORTRAN )

 ENTREES :
 --------- 
 ncf     : numero de la couleur de la face
 nca     : numero de la couleur des aretes 
 n       : nombre de sommets du polygone avec le premier=dernier sommet
 pts     : coordonnees pixels (INTEGER*2 !) des n sommets du polygone
           le pt 1 et le pt n ont memes coordonnees
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{   
  int npts, i , j;
  float counbm , courgbm[3] ;

  proc(xvcouleur)( ncf );
  XFillPolygon(display_mef, mempx, gc_mef, pts, *n, Complex, CoordModeOrigin);

  /* Traitement du postscript */
  courgbm[0] = courgb[0] ;
  courgbm[1] = courgb[1] ;
  courgbm[2] = courgb[2] ;
  counbm     = counb     ;

  proc(xvcouleur)( nca );
  XDrawLines( display_mef, mempx, gc_mef, pts, *n, CoordModeOrigin );
  /* proc(xvvoir)(); */

  /* Traitement du postscript :
     instruction ps : FP
     destination : TEMPORAIRE.EPS si l'instruction est de type dessin
                 : chaine   si l'instruction correspond a un menu */
  if (lasopsc > 0){ 
    npts = *n-1;
    iFP  = 1 ;
    if (nbrcon > 0){   
      fprintf(fpo,"%s",concat) ;
      nbrcon  = 0 ;
      concat[0] = '\0' ;
    }
    buf[0] = '\0' ;
    for (j = 0 ; j <= npts / 16 ; j++)
    {
      for (i = j * 16 ; i <= min( npts-1 , (j+1) * 16 -1 ) ; i++) 
        sprintf(&buf[strlen(buf)], "%6i %6i " , pts[i].x , ypixels-pts[i].y) ;
      if (j == npts / 16) {
        if ( strlen(buf) > 207 ) {
          sprintf(&buf[strlen(buf)],"\n") ;
          if (lasopsc < 3) {
            fprintf(fpo,"%s",buf) ;
          }
          else {
            sprintf(&chaine[lasopsc-4][strlen(chaine[lasopsc-4])],"%s",buf) ;
          }
          buf[0] = '\0' ;
        }
        if (counb != -1) {
          sprintf(&buf[strlen(buf)],"%3i %4.2f %4.2f %4.2f %4.2f ",npts,courgb[0],courgb[1],courgb[2],counb) ;
        }
        else {
          sprintf(&buf[strlen(buf)],"%3i %4.2f %4.2f %4.2f 0.00 ",npts,courgb[0],courgb[1],courgb[2]) ;
        }
        if (counbm != -1) {
          sprintf(&buf[strlen(buf)],"%4.2f %4.2f %4.2f %4.2f FP\n",courgbm[0],courgbm[1],courgbm[2],counbm) ;
        }
        else {
          sprintf(&buf[strlen(buf)],"%4.2f %4.2f %4.2f 1.00 FP\n",courgbm[0],courgbm[1],courgbm[2]) ;
        }
      }
      else {
        sprintf(&buf[strlen(buf)],"\n") ;
      }
      if (lasopsc < 3) {
        fprintf(fpo,"%s",buf) ;
      }
      else {
        sprintf(&chaine[lasopsc-4][strlen(chaine[lasopsc-4])],"%s",buf) ;
      }
      buf[0] = '\0' ;
    }
  }
} 

void proc(xvsouris)( int *notypeevent, int *nbc, int *x1, int *y1 )  
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
 BUT :  RETOURNER le type et les caracteristiques d'un EVENEMENT
 -----      ( appelable en FORTRAN )

 SORTIES :
 ---------
 notypeevent: = 0 Si ABANDON  par frappe de la touche Echappement ou @
                  supprime 27/10/2010->ou demande par clic du bouton 2 de la souris
              = 1 Si CLIC ENFONCE et RELACHE D'UN BOUTON DE LA SOURIS => x1 y1
              =-1 Si CLIC SEULEMENT ENFONCE  D'UN BOUTON DE LA SOURIS => x1 y1
              =-2 Si le pointeur de la souris a bouge                 => x1 y1
              = 2 Si FRAPPE D'UN CARACTERE AU CLAVIER  
 nbc        : seulement actif si notypeevent est non nul
              si notypeevent=-+1 nbc=numero du bouton 1 ou 2 ou 3
              si notypeevent= -2 nbc=0 (pas de bouton designe)
              si notypeevent= +2 nbc=numero du caractere dans la table ASCII
 x1,y1      : seulement actif si notypeevent=+-1 ou -2
              coordonnees pixels du point clique par rapport au coin 
              superieur gauche de la fenetre
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
 MODIF  : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE LJLL PARIS         OCTOBRE 2003
12345X7..............................................................012345678*/ 
{
  XEvent   event;
  int      flag;
  int      nb;
  char     buffer[20];
  KeySym  *keysym;   /*  incidence */
   
  *notypeevent = 0;
  flag=0;

  while (!flag)
  {
    XNextEvent(display_mef, &event);

    if( event.type == MotionNotify )
    {
      /*Deplacement de la souris
        nombre d'evenements dans la queue*/
      if( XEventsQueued( display_mef, QueuedAfterFlush ) <= 0 )
	{
	  *notypeevent = -2;         /*c'est le dernier evenement*/
	  *x1 = event.xbutton.x;
	  *y1 = event.xbutton.y; 
	  if      (event.xbutton.button == Button1) {*nbc = 1;}
	  else if (event.xbutton.button == Button2) {*nbc = 2;}
	  else if (event.xbutton.button == Button3) {*nbc = 3;}
	  else                                      {*nbc = 0;}
	  flag=1;
	}
      /*printf("xvsouris MotionNotify *notypeevent %d *nbc %d\n",*notypeevent,*nbc);*/
    }

    else if(event.type == ButtonPress)
    {
      /*Un bouton a ete enfonce, presse et non encore relache
        nombre d'evenements dans la queue*/
      if( XEventsQueued( display_mef, QueuedAfterFlush ) <= 0 )
	{
	  *notypeevent = -1;         /*c'est le dernier evenement*/
	  *x1 = event.xbutton.x;
	  *y1 = event.xbutton.y; 
	  if      (event.xbutton.button == Button1) {*nbc = 1;}
	  else if (event.xbutton.button == Button2) {*nbc = 2;}
	  else if (event.xbutton.button == Button3) {*nbc = 3;}
	  else                                      {*nbc = 0;}
	  flag=1;
	}
      /*printf("xvsouris ButtonPress *notypeevent %d *nbc %d\n",*notypeevent,*nbc);*/
    }

    else if(event.type == ButtonRelease) /*Un bouton a ete relache*/
    { 
      *notypeevent = 1;
      *x1 = event.xbutton.x;
      *y1 = event.xbutton.y; 
      if      (event.xbutton.button == Button1) {*nbc = 1;}
      else if (event.xbutton.button == Button2) {*nbc = 2;} /* 27/10/2010 *notypeevent=0;} */
      else if (event.xbutton.button == Button3) {*nbc = 3;}
      else                                      {*nbc = 0;}
      flag=1;
      /*printf("xvsouris ButtonRelease *notypeevent %d *nbc %d\n",*notypeevent,*nbc);*/
    }

    else if(event.type == KeyPress) /*Une touche du clavier a ete frappee*/
    {   
      *notypeevent = 2;
      *x1    = event.xkey.x;
      *y1    = event.xkey.y; 
      *nbc   = event.xkey.keycode ;
      keysym = 0;
      nb     = XLookupString(&event.xkey,buffer,19,keysym,NULL);
      if (nb != 0)
      {
	*nbc = buffer[0];  /*codage entier du caractere*/
	if( *nbc == 27 ) {*notypeevent=0;} /*caractere 'Echappement' => ABANDON*/
	if( *nbc == 64 ) {*notypeevent=0;} /*caractere '@'           => ABANDON*/
	flag= 1;
	/*printf("xvsouris KeyPress *notypeevent %d *nbc %d\n",*notypeevent,*nbc);*/
      }
    }
  }
}
  
void proc(xvsouris2)( int *items, int *pmin0,
                      int *notypeevent, int *ibutton, int *x1, int *y1 )  
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
 BUT :  RETOURNER le type et les caracteristiques d'un EVENEMENT
 -----  VARIANTE de xvsouris avec l'accrochage aux items traces dans la fenetre
                  ( appelable en FORTRAN )
 ENTREES:
 --------
 items : tableau des items (P ou L ou S ou V ou O) de poignees dans la fenetre
 pmin0 : -2 pour non initialise
         ou indice dans items du x pixel de la mempxaccro allume

 MODIFIE:
 --------
 pmin0 : indice dans items du x pixel de l'item le plus proche


 SORTIES :
 ---------
 notypeevent: =0 Si PAS D'EVENEMENT
              =1 Si CLIC D'UN BOUTON DE LA SOURIS   
              =2 Si FRAPPE D'UN CARACTERE AU CLAVIER 
              =5 UN BOUTON PRESSE AVEC OU SANS DEPLACEMENT
 ibutton    : seulement actif si notypeevent est non nul
              si notypeevent=1 ibutton=numero du bouton
              si notypeevent=2 ibutton=numero du caractere dans la table ASCII
 x1,y1      : seulement actif si notypeevent=1
              coordonnees pixels du point clique par rapport au coin 
              superieur gauche de la fenetre
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS              OCTOBRE 2000
12345X7..............................................................012345678*/ 
{
  XEvent  event;
  int     flag;
  int     nb,nbitem,mots,d,dmin,pmin,p,k;
  char    buffer[20];
  KeySym *keysym;   /*  incidence */
   
  *notypeevent = 0;
  flag = 0;

  while (!flag)
  {
    XNextEvent(display_mef, &event);
    if(event.type == MotionNotify || event.type == ButtonPress) 
    {
      /* Un bouton presse avec ou sans deplacement */
      *notypeevent = 5;   /* accrochage actif */
      *x1 = event.xbutton.x;
      *y1 = event.xbutton.y; 
      if (event.xbutton.button == Button1)
	*ibutton = 1;
      else if (event.xbutton.button == Button2)
	*ibutton = 2;
      else 
	*ibutton = 3;

      /* la recherche du point poignee le plus proche    */
      mots   = items[0];  /* le nombre de mots par items */
      nbitem = items[2];  /* le nombre actuel d'items    */
      dmin   = 40000;
      pmin   = -1;
      p      =  0;

      for (k=0; k<nbitem; k++)
	{
          p = p + mots;   /* x et y du pixel de l'item dans la fenetre */
	  d = (items[p]  -*x1) * (items[p]  -*x1)
            + (items[p+1]-*y1) * (items[p+1]-*y1);
	  if ( d<dmin )
	    {
	      dmin=d;
	      pmin=p;
	    }
	}

      if ( (*pmin0 != pmin)  &&  (*pmin0 >= 0) )  /* on efface l'ancien mempxaccro plus proche deja trace' */
	{
	  /* XSetFunction( display_mef, gc_mef, GXorInverted ); */
	  /* XSetFunction( display_mef, gc_mef, GXequiv );      */
	  /* XSetFunction( display_mef, gc_mef, GXxor );        */
	  XSetFunction( display_mef, gc_mef, GXorInverted );
	  XCopyArea( display_mef, mempxaccro, fenetre_mef, gc_mef, 
		     0, 0, lmempxaccro, hmempxaccro,                               /* source      */
		     items[*pmin0]-lmempxaccro/2, items[*pmin0+1]-hmempxaccro/2 ); /* destination */
	  *pmin0=-2;
	}

      if ( pmin >= 0 ) /* existence d'un item suffisamment proche */
	{

	  /* on trace le nouveau carre le plus proche */
	  XSetFunction( display_mef, gc_mef, GXand );
	  XCopyArea( display_mef, mempxaccro, fenetre_mef, gc_mef, 
		     0, 0, lmempxaccro, hmempxaccro,                             /* source */
		     items[pmin]-lmempxaccro/2, items[pmin+1]-hmempxaccro/2 );   /* destination */
          /* sauvegarde du numero d'item pour l'effacer en retour */
	  *pmin0=pmin;

	  flag=1;  /* retour demande vers l'appel de xvsouris2 */
	}
    }

    else if(event.type == ButtonRelease)
    { 
      /* l'eventuel cercle d'accrochage est efface */
      if (pmin0 != NULL)  /*if (pmin0 >= 0) corrige pour HP on efface l'ancien cercle plus proche */
	{
	  XSetFunction( display_mef, gc_mef, GXorInverted );
	  XCopyArea( display_mef, mempxaccro, fenetre_mef, gc_mef,   
		     0, 0, lmempxaccro, hmempxaccro,                               /* source */
		     items[*pmin0]-lmempxaccro/2, items[*pmin0+1]-hmempxaccro/2 ); /* destination */
	}
      /* type de l'evenement */
      *notypeevent = 1;
      *x1 = event.xbutton.x;
      *y1 = event.xbutton.y; 
      if (event.xbutton.button == Button1)
        *ibutton=1;
      else if (event.xbutton.button == Button2)
	*ibutton=2 ;
      else 
	*ibutton=3;
      flag=1;
    }

    else if(event.type == KeyPress)
    {   
      *notypeevent = 2;
      *x1 = event.xkey.x;
      *y1 = event.xkey.y; 
      *ibutton= event.xkey.keycode ;
      keysym=0;
      nb=XLookupString(&event.xkey,buffer,19,keysym,NULL);
      if (nb != 0)
      {
	*ibutton = buffer[0];
	flag= 1;
      }
    }
  }
}

void proc(deplsouris)(int *x, int *y) 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
 BUT :     Deplacer le pointeur de souris en (x,y) de la fenetre
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  XWarpPointer( display_mef, None, fenetre_mef, 0, 0, 0, 0, *x, *y );
}

void proc(xvvoir)()  
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
 BUT :     VOIR les traces encore dans la memoire tampon
 -----     ( appelable en FORTRAN )

ATTENTION: 21/1/99 
Si xvvoir est trop souvent appele cf les proc(xvvoir) de ce fichier
cela entraine un lourd trafic entre X et le window manager qui doit savoir
l'empilement des fenetres pour mettre celle de Mefisto au dessus
ce qui bloque certains wm tels kde fvwm ... mais pas certains autres 
sur ibm dec hp sun ... et toutes les fenetres sont bloquees!!!
Pour y remedier, tous les appels de proc(xvvoir) ont ete supprimes
en fait, seulement mis en commentaires
Sur hp et pc rien ne change en execution si ce n'est que c'est plus rapide
et SURTOUT sans blocage avec kde fvwm ... sur PC
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  XRaiseWindow (display_mef,fenetre_mef);
  XFlush(display_mef);
}


void proc(xvpause)()  
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :  FAIRE une PAUSE jusqu'a l'entree d'un caractere au clavier 
 -----    ( appelable en FORTRAN )
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  XEvent event;

  XNextEvent(display_mef, &event);
  switch(event.type)
  {
    case KeyPress: break;
  }
}

                                                
void proc(xvfbordrectangle)( int *x, int *y, int *width, int *height )  
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :     TRACER les aretes du contour d'un rectangle dans fenetre_mef
 -----       ( appelable en FORTRAN )

 ENTREES :
 ---------  
 x,y     : coordonnees pixels du coin superieur gauche
 width   : largeur en pixels du rectangle
 height  : hauteur en pixels du rectangle
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  XDrawRectangle(display_mef, fenetre_mef, gc_mef, *x, *y, *width, *height);
}

void proc(xvbordrectangle)( int *x, int *y, int *width, int *height )  
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :     TRACER les aretes du contour d'un rectangle
 -----       ( appelable en FORTRAN )

 ENTREES :
 ---------  
 x,y     : coordonnees pixels du coin superieur gauche
 width   : largeur en pixels du rectangle
 height  : hauteur en pixels du rectangle
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  XDrawRectangle(display_mef, mempx, gc_mef, *x, *y, *width, *height);
 /*  proc(xvvoir)(); */

  /* Traitement du postscript :
     instruction ps : r
     destination : TEMPORAIRE.EPS si l'instruction est de type dessin
                 : chaine   si l'instruction correspond a un menu */
  if (lasopsc > 0){ 
    ire = 1 ;
    if (nbrcon > 0){   
      fprintf(fpo,"%s",concat) ;
      nbrcon  = 0 ;
      concat[0] = '\0' ;
    }
    buf[0] = '\0' ;
    if (counb != -1) {
      sprintf(&buf[strlen(buf)], "%6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f r\n",
              *width, -*height, *x, ypixels-*y,courgb[0],courgb[1],courgb[2], counb) ;
    }
    else {
      sprintf(&buf[strlen(buf)], "%6i %6i %6i %6i %4.2f %4.2f %4.2f 0.00 r\n",
              *width, -*height, *x, ypixels-*y,courgb[0],courgb[1],courgb[2]) ;
    }
    if (lasopsc < 3){
      fprintf(fpo,"%s",buf) ;
    }
    else{
      sprintf(&chaine[lasopsc-4][strlen(chaine[lasopsc-4])],"%s",buf) ;
    }
  }
} 

void proc(xvfrectangle)(int *x,int *y,int *width,int *height)  
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :     TRACER le remplissage du contour d'un rectangle dans la fenetre_mef
 -----       ( appelable en FORTRAN )

 ENTREES :
 ---------  
 x,y     : coordonnees pixels du coin superieur gauche
 width   : largeur en pixels du rectangle
 height  : hauteur en pixels du rectangle   vers le bas
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  XFillRectangle(display_mef, fenetre_mef, gc_mef, *x, *y, *width, *height); 
}
 

void proc(xvrectangle)(int *x,int *y,int *width,int *height)  
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :     TRACER le remplissage du contour d'un rectangle
 -----       ( appelable en FORTRAN )

 ENTREES :
 ---------  
 x,y     : coordonnees pixels du coin superieur gauche
 width   : largeur en pixels du rectangle
 height  : hauteur en pixels du rectangle   vers le bas
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  XFillRectangle(display_mef, mempx, gc_mef, *x, *y, *width, *height); 
 /*  proc(xvvoir)(); */

  /* Traitement du postscript :
     instruction ps : R
     destination : TEMPORAIRE.EPS si l'instruction est de type dessin
                 : chaine   si l'instruction correspond a un menu */
  if (lasopsc > 0){              
    iRe = 1 ;
    if (nbrcon > 0){   
      fprintf(fpo,"%s",concat) ;
      nbrcon  = 0 ;
      concat[0] = '\0' ;
    }
    buf[0] = '\0' ;
    if (counb != -1) {
      sprintf(&buf[strlen(buf)], "%6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f R\n",
              *width, -*height, *x, ypixels-*y,courgb[0],courgb[1],courgb[2],counb) ;
    }
    else {
      sprintf(&buf[strlen(buf)], "%6i %6i %6i %6i %4.2f %4.2f %4.2f 1.00 R\n",
              *width, -*height, *x, ypixels-*y,courgb[0],courgb[1],courgb[2]) ;
    }
    if (lasopsc < 3){
      fprintf(fpo,"%s",buf) ;
    }
    else{
      sprintf(&chaine[lasopsc-4][strlen(chaine[lasopsc-4])],"%s",buf) ;
    }
  }
}
  
   
void proc(xvbordarcellipse)( int *x, int *y, int *width, int *height,
                             float *angle1, float *angle2 )  
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :     TRACER le bord d'un secteur d'une ellipse
 -----       ( appelable en FORTRAN )

 ENTREES :
 ---------  
 x,y     : coordonnees pixels du centre de l'ellipse
 width   : demi-largeur en pixels de l'ellipse
 height  : demi-hauteur en pixels de l'ellipse
 angle1  : angle de depart (a partir de l'axe Ox sens direct en degres )
 angle2  : angle du secteur a tracer (a partir de l'angle de depart en degres )
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{ 
  int adep , afin ;
  adep = (int) (*angle1 * 64) ; 
  afin = (int) (*angle2 * 64) ; 
  XDrawArc(display_mef, mempx, gc_mef, *x-*width, *y-*height, *width * 2,
           *height * 2, adep , afin );
 /*  proc(xvvoir)(); */

  /* Traitement du postscript :
     instruction ps : el
     destination : TEMPORAIRE.EPS si l'instruction est de type dessin
                 : chaine   si l'instruction correspond a un menu */
  if (lasopsc > 0){              
    iel = 1 ;
    if (nbrcon > 0){   
      fprintf(fpo,"%s",concat) ;
      nbrcon  = 0 ;
      concat[0] = '\0' ;
    }
    buf[0] = '\0' ;
    if (*angle2 >= 0){
      adep = *angle1 ;  
      afin = *angle1+*angle2 ;  
    }
    else{
      afin = *angle1 ;  
      adep = *angle1+*angle2 ;  
    }
    if (counb != -1) {
      sprintf(&buf[strlen(buf)], "%6i %6i %6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f el\n",
              adep , afin , *width, *height, *x, ypixels-*y,courgb[0],courgb[1],courgb[2], counb) ;
    }
    else {
      sprintf(&buf[strlen(buf)], "%6i %6i %6i %6i %6i %6i %4.2f %4.2f %4.2f 0.00 el\n",
              adep , afin , *width, *height, *x, ypixels-*y,courgb[0],courgb[1],courgb[2]) ;
    }
    if (lasopsc < 3){
      fprintf(fpo,"%s",buf) ;
    }
    else{
      sprintf(&chaine[lasopsc-4][strlen(chaine[lasopsc-4])],"%s",buf) ;
    }
  }
}
   

void proc(xvarcellipse)( int *x, int *y, int *width, int *height,
                         float *angle1, float *angle2 ) 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :     TRACER le remplissage d'un secteur d'une ellipse
 -----       ( appelable en FORTRAN )

 ENTREES :
 ---------  
 x,y     : coordonnees pixels du centre de l'ellipse
 width   : demi-largeur en pixels de l'ellipse
 height  : demi-hauteur en pixels de l'ellipse
 angle1  : angle de depart (a partir de l'axe Ox sens direct en degres )
 angle2  : angle du secteur a tracer (a partir de l'angle de depart en degres )
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS                  MAI 1994
12345X7..............................................................012345678*/ 
{
  int adep , afin ;
  adep = (int) (*angle1 * 64) ; 
  afin = (int) (*angle2 * 64) ; 
  XFillArc(display_mef, mempx, gc_mef, *x-*width, *y-*height, *width * 2,
           *height * 2, adep , afin );
 /*  proc(xvvoir)(); */

  /* Traitement du postscript :
     instruction ps : El
     destination : TEMPORAIRE.EPS si l'instruction est de type dessin
                 : chaine   si l'instruction correspond a un menu */
  if (lasopsc > 0){              
    iEl = 1 ;
    if (nbrcon > 0){   
      fprintf(fpo,"%s",concat) ;
      nbrcon  = 0 ;
      concat[0] = '\0' ;
    }
    buf[0] = '\0' ;
    if (*angle2 >= 0){
      adep = *angle1 ;  
      afin = *angle1+*angle2 ;  
    }
    else{
      afin = *angle1 ;  
      adep = *angle1+*angle2 ;  
    }
    if (counb != -1) {
      sprintf(&buf[strlen(buf)], "%6i %6i %6i %6i %6i %6i %4.2f %4.2f %4.2f %4.2f El\n",
              adep , afin , *width, *height, *x, ypixels-*y,courgb[0],courgb[1],courgb[2], counb) ;
    }
    else {
      sprintf(&buf[strlen(buf)], "%6i %6i %6i %6i %6i %6i %4.2f %4.2f %4.2f 1.00 El\n",
              adep , afin , *width, *height, *x, ypixels-*y,courgb[0],courgb[1],courgb[2]) ;
    }
    if (lasopsc < 3){
      fprintf(fpo,"%s",buf) ;
    }
    else{
      sprintf(&chaine[lasopsc-4][strlen(chaine[lasopsc-4])],"%s",buf) ;
    }
  }
}        


void proc(tempscpu) ( double *tclock )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :   Retourne le temps CPU utilise 
 -----
 SORTIES :
 --------- 
 tclock  : le temps CPU utilise
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AUTEUR : PERRONNET ALAIN    ANALYSE NUMERIQUE UPMC  PARIS          NOVEMBRE 1994
12345X7..............................................................012345678*/  
{  
  *tclock = ( (double) clock() ) / ( (double) CLOCKS_PER_SEC );
}



void proc(secondes1970) ( double *secondes )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :   Retourne le nombre de secondes depuis le 1/1/1970 minuit
 -----
 SORTIES :
 --------- 
 secondes : le nombre de secondes depuis le 1/1/1970 minuit
            AVEC les micro-secondes factices ont pour but de dissocier 2 appels
            memes tres rapproches et notamment dans la meme seconde
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AUTEUR : PERRONNET ALAIN    ANALYSE NUMERIQUE UPMC  PARIS              MARS 1996
12345X7..............................................................012345678*/  
{  
   nbs1970 = nbs1970 + 1 ;
   *secondes = (double) time( NULL ) + ( (double) nbs1970 ) * 0.000001 ;
}  


void proc(secondes1969) ( double *secondes )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :   Retourne le nombre de secondes depuis le 1/1/1970 minuit
 -----
 SORTIES :
 --------- 
 secondes : le nombre de secondes depuis le 1/1/1970 minuit
            SANS les micro-secondes factices de secondes1970
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AUTEUR : ALAIN PERRONNET LJLL UPMC & ST PIERRE DU PERRAY            OCTOBRE 2010
12345X7..............................................................012345678*/  
{  
   *secondes = (double) time( NULL ) ;
}


void proc(nomordinateurhote) ( char *host, int *nbcar)
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :   Retourne le nom de l'ordinateur hote
 -----
 SORTIES :
 --------- 
 host    : nom de l'ordinateur hote
 nbcar   : nombre de caracteres du nom de l'ordinateur hote
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AUTEUR : P LAUG ou F. HECHT INRIA                                  DECEMBRE 1994
12345X7..............................................................012345678*/  
{
   char name[80];

   if (gethostname(name, 80) != 0)
   {
      fprintf(stderr, "nom de l'ordinateur hote a declarer\n") ;
      return ;
   }
   *nbcar = strlen(name) ;
   strcpy( host, name );
}



void proc(ladate) ( int *a, int *m, int *j )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :   Retourne la date du jour en an mois jour
 -----
 SORTIES :
 --------- 
 a       : annee 
 m       : mois
 j       : jour
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AUTEUR : P LAUG ou F. HECHT INRIA                                  DECEMBRE 1994
12345X7..............................................................012345678*/  
{
    time_t timevar;
    struct tm *pttm;

    timevar=time(&timevar);
    pttm = localtime(&timevar);
 
    *a = pttm->tm_year;
    *m = pttm->tm_mon+1;
    *j = pttm->tm_mday ;
}



void proc(heureminuteseconde) ( int *h, int *m, int *s, int *millis )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :   Retourne l'heure en heure minute seconde et millisecondes
 -----
 SORTIES :
 --------- 
 h       : heure
 m       : minute
 s       : seconde
 millis  : millisecondes
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AUTEUR : P LAUG ou F. HECHT INRIA                                  DECEMBRE 1994
12345X7..............................................................012345678*/  
{
    time_t timevar;
    struct tm *pttm;

    timevar=time(&timevar);
    pttm = localtime(&timevar);
    *h = pttm->tm_hour;
    *m = pttm->tm_min;
    *s = pttm->tm_sec;
    *millis = 0;
}
      

void proc(valvarenv)( char *nom, int *lval_admis,
                      char *val, int *lval_trouve )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 BUT :   Retourne la valeur de la variable d'environnement nom pour le Fortran
 -----   Attention a terminer le nom par CHAR(0)

         Exemple d'appel :
         CHARACTER*10  C
         INTEGER       L
         CALL VALVARENV( 'DISPLAY'//CHAR(0), LEN(C),  C, L )

 ENTREES :
 ---------
 nom        : nom de la variable d'environnement
 lval_admis : nombre maximal de caracteres admis pour sa valeur

 SORTIES :
 --------- 
 val         : la chaine de caracteres valeur de la variable d'environnement
 lval_trouve : < 0      si le nom de la variable n'est pas trouvee
               > LEN(C) si la valeur de la variable a ete tronquee
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AUTEUR : P LAUG ou F. HECHT INRIA                                      Mars 1991
12345X7..............................................................012345678*/  
{
   char  *ptenv, *getenv() ;
   int    i ;

   ptenv = getenv(nom) ;
   if (ptenv == NULL) 
   {
      *lval_trouve = -1 ;  /* variable d'environnement non trouvee */
      return ;
   }
   for (i=0 ; *(ptenv+i)!='\0' ; i++)
      if (i<*lval_admis) *(val+i) = *(ptenv+i) ;
   *lval_trouve = i ;
}

void proc(xvinitierps)( int *modeps )
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :   INITIALISER le POSTSCRIPT
 -----       ( appelable en FORTRAN )
 ENTREE  :
 ---------  
 modeps : 0 mode standard bibliotheque xvue
          1 mode mefisto
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : DOURSAT CHRISTOPHE UPMC ANALYSE NUMERIQUE PARIS              JUIN 1994
12345X7..............................................................012345678*/ 
{
  int lasops ;
  lasopsc = 0;
  modepsc = *modeps ;
  icolorm = ndcoul + 1 ;
  counb      = 0.0 ;
  courgb[0] = 0.0 ;  courgb[1] = 0.0 ;  courgb[2] = 0.0 ;
  /* initialisation des pointeurs et variables postscript */
  fpi = NULL ; 
  fpo = NULL ; 
  menu = 0 ;
  if ( modepsc != 0 ){
  /* Initialisation des tailles des chaines-tampons des menus
     Nombre de caracteres a reserver par instruction PS
       1 rectangle                  = 35
       1 segment                    = 35
       1 changement epaisseur trait = 10
       1 chaine texte               = 24 + long chaine
     Par menu il y a :
       2     rectangles (un fond et un contour)
       2     segments de bordure
       3     changements d'epaisseur (maximum)
       n     chaines de caracteres
       (n-1) segments entre les chaines de caracteres 
       + 50  de securite */

  /* Documentation */
  longchaine[0] = 2*35 + 2*35 + 3*10 + MXLGDO*(NBCADO+24) + (MXLGDO-1) * 35 + 50 ;   
  /* Historique */
  longchaine[1] = 2*35 + 2*35 + 3*10 + MXLGHI*(NBCAHI+24) + (MXLGHI-1) * 35 + 50 ;
  /* Menu */
  longchaine[2] = 2*35 + 2*35 + 3*10 + MXLGME*(NBCAME+24) + (MXLGME-1) * 35 + 50 ;
  /* Lignes lues */
  longchaine[3] = 2*35 + 2*35 + 3*10 + MXKLG *(NBCALI+24)                   + 50 ;
  /* Invite */
  longchaine[4] = 2*35 + 2*35 + 3*10 + MXLGIN*(NBCAIN+24) + (MXLGIN-1) * 35 + 50 ;
  /* Ligne de saisie */
  longchaine[5] = 2*35 + 2*35 + 3*10 + MXLGSA*(NBCALG+24) + (MXLGSA-1) * 35 + 50 ;
  /* Erreur */
  longchaine[6] = 2*35 + 2*35 + 3*10 + MXLGER*(NBCAER+24) + (MXLGER-1) * 35 + 50 ;
  /* Histogrammes et Qualites : 1+32 rectangles et chaines de caracters de 72 */
  longchaine[7] =   35 +32*35        +     32*(  72  +24)                   + 50 ;
  }
  /* effacement si necessaire du fichier ps */
  lasops = 0 ;
  proc(xvpostscript)(&lasops) ;
  /* ouverture du futur fichier ps */
  lasops = 1 ;
  proc(xvpostscript)(&lasops) ;
}

void proc(xvimprimerps)( char nomfichier[], int *length ) 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :   IMPRIMER nomfichier.eps ( appelable en FORTRAN )
 -----

 ENTREE  :
 ---------  
 nomfichier : CHAINE DE CARACTERE
 length     : LONGUEUR UTILE DE LA CHAINE DE CARACTERE

 SORTIE  :
 ---------  
 length     : 0 si pas d'erreur
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : DOURSAT CHRISTOPHE UPMC ANALYSE NUMERIQUE PARIS              JUIN 1994
12345X7..............................................................012345678*/
{ /*int i;*/
  int nbc ;
        
  buf[0] = '\0' ;
  nbc    = *length ;
  sprintf(format,"%%.%ds.eps", nbc) ;
  sprintf(buf,format,nomfichier) ;
  /* changement du titre si possible */
/*  if ((fpi = fopen(buf,"r+"))!=NULL) {
    rewind(fpi) ;
    for (i=1 ; i<9 ; i++) { fgets(buf,255,fpi) ; }
    buf[0] = '\0' ;
    sprintf(format,"%%.%ds.eps) /titre exch def\n",nbc) ;
    sprintf(buf,format,nomfichier) ;
    fprintf(fpi,"(%s",buf) ;
    fclose(fpi) ;
  } */

/* ******************************** ATTENTION ****************************************
                  ORDRE D'IMPRESSION DEPENDANT TYPE IMPRIMANTE ET SYSTEME             
   ******************************** ATTENTION **************************************** */
   sprintf(format,"lpr -Ppastis %%.%ds.eps",nbc); /* ATTENTION ici imprimante de LJLL */

/* sprintf(format,"lp -dlaserhp  %%.%ds.eps",nbc) ; */   /* ATTENTION ici HP */
/* sprintf(format,"lpr  %%.%ds.eps",nbc) ; */  /* ajouter -Plaser ? */
/* sprintf(format,"/bin/prf -trans %%.%ds.eps",nbc) ; */  /* Apollo LaserWriterII */

/* ******************************** ATTENTION **************************************** */
   sprintf(buf,format,nomfichier) ;            
   nbc = system(buf) ;
}
       
void proc(xvsauverps)( char nomfichier[], int *length ) 
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 BUT :   CONCATENER LES DONNEES DE DESSIN PS CONTENUES DANS TEMPORAIRE.EPS
 -----   AVEC LES PROCEDURES PS ET LES DIVERS MENUS EVENTUELS
         DANS nomfichier.eps ( appelable en FORTRAN )

 ENTREE  :
 ---------  
 nomfichier : CHAINE DE CARACTERES
 length     : LONGUEUR UTILE DE LA CHAINE DE CARACTERE
 ATTENTION  : length ne doit pas etre l'ADRESSE d'UNE CONSTANTE
              car son ecrasement dans cette fonction produit SEGMENTATION FAULT!
              en FORTRAN, il faut utiliser une VARIABLE ENTIERE, PAS une CONSTANTE!

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 AUTEUR : DOURSAT CHRISTOPHE UPMC ANALYSE NUMERIQUE PARIS              JUIN 1994
12345X7..............................................................012345678*/
{ int  i;
  char cpwd[100];
  char *ligne, *car ;
  int  an,mo,jo;

  /* verification des fichiers */
  fprintf(fpo,"%s",concat) ;
  nbrcon  = 0 ;
  concat[0] = '\0' ;
  i = 0 ;
  if (fpi != NULL) {        
    fclose(fpi) ;
    fpi = NULL ; 
  } 
  sprintf(format,"%%.%ds.eps",*length) ;
  sprintf(buf,format,nomfichier) ;
  if ((fpi = fopen(buf,"w"))==NULL) {
    printf ("Erreur de creation de fichier\n") ; 
    i = 1 ;
    *length = -2 ;
  }
  if ( i == 0 ){            
    if (fpo != NULL) {        
      fclose(fpo) ;
    }                                             
    fpo = fpi  ;
    i = 0 ;
    fpi = NULL ;
    if ((fpi = fopen("TEMPORAIRE.EPS","r"))==NULL) {
      printf ("Erreur d'ouverture de fichier\n") ; 
      i = 1 ;
      *length = -2 ;
    }
  }
  if ( i == 0 ){            
    lasopsc = 0 ;
    /* ecriture du dictionnaire de procedures PS utilisees */
    fprintf(fpo,"%%!PS-Adobe\n") ;
    if ( modepsc == 0 ){
      fprintf(fpo,"%%%%Creator: xvue\n") ;
    } else { 
      fprintf(fpo,"%%%%Creator: Mefisto\n") ;
      sprintf(buf,"LOGNAME") ;
      i = 100 ;
      proc(valvarenv)( buf, &i, cpwd, &an) ;
      sprintf(format,"%%%%%%%%User: %%.%ds\n",an) ;
      for (i=0 ; (i<an) && (*(cpwd+i)!='.') ; i++) ;
      sprintf(format,"%%%%%%%%User: %%.%ds\n",i) ;
      fprintf(fpo,format,cpwd) ;
   }
    sprintf(format,"%%%%%%%%Title: %%.%ds",*length) ;
    sprintf(buf,format,nomfichier) ;
    fprintf(fpo,"%s\n",buf) ;
    proc(ladate)(&an,&mo,&jo) ;
    an = an % 100;
    fprintf(fpo,"%%%%CreationDate: %2i/%2i/%2i\n",jo,mo,an) ;
    fprintf(fpo,"%%%%BoundingBox:     0     0 %5i %5i\n",xpixels,ypixels) ;
    fprintf(fpo,"%%%%EndComments\n\n") ;
    fprintf(fpo,"%%=======================================================================================\n") ;
    fprintf(fpo,"%% A NE PAS CHANGER\n") ;
    fprintf(fpo,"%%=======================================================================================\n") ;
    fprintf(fpo,"  mark\n") ;
    fprintf(fpo,"  /tex where { pop /texval true def }{ /texval false def } ifelse\n") ;
    fprintf(fpo,"  /texdvips where { pop /texdvipsval true def /texval true def }\n") ;
    fprintf(fpo,"                  { /texdvipsval false def } ifelse\n") ;
    fprintf(fpo,"  /ptcm  28.3464566929134 def /cm {ptcm mul} def\n") ;
    fprintf(fpo,"%%=======================================================================================\n") ;
    fprintf(fpo,"%% OPTIONS UTILISATEURS\n") ;
    fprintf(fpo,"%%=======================================================================================\n") ;
    fprintf(fpo,"%% OPTION COULEUR\n") ;
    fprintf(fpo,"%%=======================================================================================\n") ;
    fprintf(fpo,"                            %% trace couleur (true) ou n&b = palette gris (false)\n") ;
    fprintf(fpo,"true   /couleur  exch def   %% remis systematiquement a true si le fichier est envoye sur\n") ;
    fprintf(fpo,"                            %% une imprimante couleur referencee (cf. Options Imprim.)\n\n") ;
    fprintf(fpo,"                            %% couleur de fond : true = noir false = blanc\n") ;
    fprintf(fpo,"false  /fondcou  exch def   %% remis systematiquement a false si couleur = false\n") ;
    fprintf(fpo,"                            %% sauf si imprimante couleur referencee\n\n") ;
    fprintf(fpo,"%% RECADRAGE\n") ;
    fprintf(fpo,"%%=======================================================================================\n") ;
    fprintf(fpo,"    0  /zoomxmin exch def   %% \n") ;
    fprintf(fpo,"%5i  /zoomxmax exch def   %% definition du rectangle de zoom choisi\n",xpixels) ;
    fprintf(fpo,"    0  /zoomymin exch def   %% a utiliser pour recadrer la figure\n") ;
    fprintf(fpo,"%5i  /zoomymax exch def   %%\n\n",ypixels) ;
    fprintf(fpo,"                            %% encadrement automatique des zones de trace\n") ;
    fprintf(fpo,"false  /autoreca exch def   %% taper dans le fichier autoreca pour definir la limite\n") ;
    fprintf(fpo,"                            %% entre deux zones a encadrer\n\n") ;
    fprintf(fpo,"%% OPTIONS DE TRACE\n") ;
    fprintf(fpo,"%%=======================================================================================\n") ;
    fprintf(fpo,"                            %% largeur du trait de trace ]0,+oo[, 1 etant la norme.\n") ;
    fprintf(fpo,"   -1  /lartrait exch def   %% le code -1 (valeur par defaut) permet de conserver la meme\n") ;
    fprintf(fpo,"                            %%   largeur de trait et d'ecriture quelque soit le zoom.\n\n") ;
    if ( modepsc != 0 ){
      fprintf(fpo,"                            %% option de trace que maillage         : 0\n") ;
      if ( menu ){
        fprintf(fpo,"    2  /tracemen exch def   %%                 avec legende qualite : 1\n") ;
      } else {
        if ( strlen(chaine[7]) != 0 ){
          fprintf(fpo,"    1  /tracemen exch def   %%                 avec legende qualite : 1\n") ;
        } else {
          fprintf(fpo,"    0  /tracemen exch def   %%                 avec legende qualite : 1\n") ;
        }
      }
      fprintf(fpo,"                            %%                 avec menus           : 2\n") ;
      fprintf(fpo,"                            %% par defaut : si il y a un menu 2\n") ;
      fprintf(fpo,"                            %%              sinon : si qualite 1, 0 sinon\n\n") ;
    }
    fprintf(fpo,"                            %% option de trace des facettes 0 fil de fer\n") ;
    fprintf(fpo,"    2  /tracecou exch def   %%                              1 faces blanches\n") ;
    fprintf(fpo,"                            %%                              2 faces de couleur\n\n") ;
    fprintf(fpo," true  /traceecr exch def   %% option de trace des chaines de caracteres\n\n") ;
    fprintf(fpo,"%% OPTIONS DE MISE EN PAGE\n") ;
    fprintf(fpo,"%%=======================================================================================\n") ;
    fprintf(fpo," true  /tracecad exch def   %% trace d'un cadre true ou false\n") ;
    fprintf(fpo,"                            %% \n") ;
    fprintf(fpo,"    4  /epcadre  exch def   %% epaisseur du cadre\n\n") ;
    if ( modepsc != 0 ){
      if ( menu ){
        fprintf(fpo,"false  /tracelog exch def   %% trace du logo Mefisto : par defaut false si menu true sinon\n") ;
      } else {
        fprintf(fpo," true  /tracelog exch def   %% trace du logo Mefisto : par defaut false si menu true sinon\n") ;
      }
    }
    fprintf(fpo,"    1  /xpolog   exch def   %% position horizontale du logo Mefisto -1 gauche 1 droite\n") ;
    fprintf(fpo,"   -1  /ypolog   exch def   %% position verticale   du logo Mefisto -1 bas    1 haut\n") ;
    fprintf(fpo,"                            %% position par defaut 1 -1 coin bas droit\n\n") ;

/*  *********************************** ATTENTION ***************************************************************
                          DONNEES DEPENDANTES TYPE IMPRIMANTE
    *********************************** ATTENTION ************************************************************ */

    fprintf(fpo,"%% OPTIONS DEPENDANTES DE L'IMPRIMANTE\n") ;
    fprintf(fpo,"%%=======================================================================================\n") ;
    fprintf(fpo,"texval not {\n") ;
    fprintf(fpo,"  27.0 cm  /vsize exch def  %% hauteur max de la zone de trace sur A4 < 29.7cm\n") ;
    fprintf(fpo,"  19.0 cm  /hsize exch def  %% largeur max de la zone de trace sur A4 < 21.0cm\n") ;
    fprintf(fpo,"} if\n\n") ;
    fprintf(fpo,"version (2010.113) eq { %% palette grise pour LaserWriter II (salle apollo)\n") ;
    fprintf(fpo,"  0.0  /gsefo    exch def   %% niveau du gris le plus fonce pour segment 0 = noir\n") ;
    fprintf(fpo,"  0.8  /gsecl    exch def   %%                        clair              1 = blanc\n") ;
    fprintf(fpo,"  0.1  /gfafo    exch def   %% niveau du gris le plus fonce pour faces   0 = noir\n") ;
    fprintf(fpo,"  1.0  /gfacl    exch def   %%                        clair              1 = blanc\n") ;
    fprintf(fpo,"}{\n") ;
    fprintf(fpo,"version (2011.110) eq { %% palette grise pour LaserJet4 (salle HP)\n") ;
    fprintf(fpo,"  0.1  /gsefo    exch def   %% niveau du gris le plus fonce pour segment 0 = noir\n") ;
    fprintf(fpo,"  0.9  /gsecl    exch def   %%                        clair              1 = blanc\n") ;
    fprintf(fpo,"  0.5  /gfafo    exch def   %% niveau du gris le plus fonce pour faces   0 = noir\n") ;
    fprintf(fpo,"  1.0  /gfacl    exch def   %%                        clair              1 = blanc\n") ;
    fprintf(fpo,"}{\n") ;
    fprintf(fpo,"version (xxxx.xxx) eq { %% palette grise pour LaserWriterII (salle Mac)\n") ;
    fprintf(fpo,"  0.1  /gsefo    exch def   %% niveau du gris le plus fonce pour segment 0 = noir\n") ;
    fprintf(fpo,"  0.9  /gsecl    exch def   %%                        clair              1 = blanc\n") ;
    fprintf(fpo,"  0.5  /gfafo    exch def   %% niveau du gris le plus fonce pour faces   0 = noir\n") ;
    fprintf(fpo,"  1.0  /gfacl    exch def   %%                        clair              1 = blanc\n") ;
    fprintf(fpo,"}{\n") ;
    fprintf(fpo,"version (2013.104) eq { %% couleur pour jet d'encre couleur\n") ;
    fprintf(fpo,"  /couleur true def\n") ;
    fprintf(fpo,"}{                      %% palette grise par defaut\n") ;
    fprintf(fpo,"  0.0  /gsefo    exch def   %% niveau du gris le plus fonce pour segment 0 = noir\n") ;
    fprintf(fpo,"  0.9  /gsecl    exch def   %%                        clair              1 = blanc\n") ;
    fprintf(fpo,"  0.5  /gfafo    exch def   %% niveau du gris le plus fonce pour faces   0 = noir\n") ;
    fprintf(fpo,"  1.0  /gfacl    exch def   %%                        clair              1 = blanc\n") ;
    fprintf(fpo,"} ifelse } ifelse } ifelse } ifelse\n\n") ;
/* ******************************** ATTENTION **************************************** */
    fprintf(fpo,"%% OPTION FONTES\n") ;
    fprintf(fpo,"%%=======================================================================================\n") ;
    fprintf(fpo,"                               %% fontes par defaut (a mettre sous la forme par ex\n") ;
    fprintf(fpo,"/Courier      /fontm0 exch def %%   de /Courier-BoldItalic /Helvetica ...)\n") ;
    fprintf(fpo,"/Helvetica    /fontp0 exch def %% * /fontm0 est la fonte monospaced\n") ;
    fprintf(fpo,"                               %% * /fontp0 est la fonte proportional\n\n") ;
    fprintf(fpo,"                               %%   de /Courier-BoldItalic /Helvetica ...)\n\n") ;
    fprintf(fpo,"false         /bavard exch def %% bavard a true indique les chaines de caracteres\n") ;
    fprintf(fpo,"                               %%   ou la police choisie n'existe pas\n\n") ;
    fprintf(fpo,"%%=======================================================================================\n\n") ;
    fprintf(fpo,"%%%%BeginProcSet \n") ;
    fprintf(fpo,"%%%%\n") ;
    fprintf(fpo,"%%%% Auto-encadrement\n") ;
    fprintf(fpo,"autoreca{ %% definition des macros servant a l'auto-encadrement\n") ;
    fprintf(fpo,"  /nbrcad 0 def\n") ;
    fprintf(fpo,"  /c@rli {currentlinewidth} def\n") ;
    fprintf(fpo,"  /initcad{/nbrcad nbrcad 1 add def /xmin{dup}def /ymin{dup}def /xmax{dup}def /ymax{dup}def}def\n") ;
    fprintf(fpo,"  initcad\n") ;
    fprintf(fpo,"  /t@stxy{2 copy 2 copy /n c@rli 2 div def\n") ;
    fprintf(fpo,"    n sub dup ymin le{/ymin exch def}{pop}ifelse n sub dup xmin le{/xmin exch def}{pop}ifelse\n") ;
    fprintf(fpo,"    n add dup ymax ge{/ymax exch def}{pop}ifelse n add dup xmax ge{/xmax exch def}{pop}ifelse}def\n") ;
    fprintf(fpo,"  /t@stx@{2 copy dup ymin le{/ymin exch def}{pop}ifelse dup xmin le{/xmin exch def}{pop}ifelse\n") ;
    fprintf(fpo,"                 dup ymax ge{/ymax exch def}{pop}ifelse dup xmax ge{/xmax exch def}{pop}ifelse}def\n") ;
    fprintf(fpo,"  /m@veto {t@stxy moveto}def\n") ;
    fprintf(fpo,"  /m@v@to {2 copy t@stx@ moveto}def\n") ;
    fprintf(fpo,"  /m@v@t@ {2 copy 2 copy 2 copy tfont add t@stx@ t@stx@ moveto}def\n") ;
    fprintf(fpo,"  /sh@w   {3 2 roll dup show stringwidth 3 2 roll add 3 1 roll add exch t@stx@}def\n") ;
    fprintf(fpo,"  /l@neto {t@stxy lineto}def\n") ;
    fprintf(fpo,"  /l@n@to {2 copy t@stx@ lineto}def\n") ;
    fprintf(fpo,"  /rl@neto{2 copy currentpoint 3 2 roll add 3 1 roll add exch t@stxy pop pop rlineto}def\n") ;
    fprintf(fpo,"  /rl@n@to{2 copy currentpoint 3 2 roll add 3 1 roll add exch t@stx@ rlineto}def\n") ;
    fprintf(fpo,"  /t@stel{6 copy /ycent exch def /xcent exch def /hautel exch def /largel exch def\n") ;
    fprintf(fpo,"    2 copy exch dup truncate 0 ge{90 idiv 1 add}{90 idiv}ifelse\n") ;
    fprintf(fpo,"    exch dup truncate 0 ge{90 idiv}{90 idiv 1 sub}ifelse\n") ;
    fprintf(fpo,"    2 copy exch sub 1 add  3 1 roll 1 exch {90 mul exch} for\n") ;
    fprintf(fpo,"    2 add {dup cos largel mul xcent add exch sin hautel mul ycent add t@stxy pop pop} repeat}def\n") ;
    fprintf(fpo,"  /t@st@l{6 copy /ycent exch def /xcent exch def /hautel exch def /largel exch def\n") ;
    fprintf(fpo,"    2 copy exch dup truncate 0 ge{90 idiv 1 add}{90 idiv}ifelse\n") ;
    fprintf(fpo,"    exch dup truncate 0 ge{90 idiv}{90 idiv 1 sub}ifelse\n") ;
    fprintf(fpo,"    2 copy exch sub 1 add  3 1 roll 1 exch {90 mul exch} for\n") ;
    fprintf(fpo,"    2 add {dup cos largel mul xcent add exch sin hautel mul ycent add t@stx@} repeat}def\n") ;
    fprintf(fpo,"  /autoreca{xmin xmax ymin ymax initcad}def\n") ;
    fprintf(fpo,"  /encadre {xmin xmax ymin ymax 1 epais\n") ;
    fprintf(fpo,"    /font0 /Courier-Bold findfont [15 0 0 17 0 0] makefont definefont pop\n") ;
    fprintf(fpo,"    lartrait -1 eq {/font0 /font0 findfont reduc scalefont definefont pop} if\n") ;
    fprintf(fpo,"    /taillefonte{ dup /tfont exch 17 div 14 mul reduc mul def 14.0 div /font0 findfont exch scalefont setfont } bind def\n") ;
    fprintf(fpo,"    17 taillefonte /hsi (1234567890) stringwidth pop def\n") ;
    fprintf(fpo,"    /tab {hsi 0 rmoveto dup stringwidth pop neg 0 rmoveto show} def\n") ;
    fprintf(fpo,"    /hsize hsi 4 mul (1 : ) stringwidth pop add def\n") ;
    fprintf(fpo,"    /vsize tfont nbrcad 1 add mul 1.2 mul def /angle 0 def\n") ;
    fprintf(fpo,"    texval {\n") ;
    fprintf(fpo,"      /xorig zoomxmin tracecad{epcadre sub}if def /yorig zoomymin tracecad{epcadre sub}if def\n") ;
    fprintf(fpo,"      /cadrepath{newpath\n") ;
    fprintf(fpo,"      xorig yorig vsize sub moveto zoomxmax tracecad{epcadre add}if yorig vsize sub lineto\n") ;
    fprintf(fpo,"      zoomxmax tracecad{epcadre add}if zoomymax tracecad{epcadre add}if lineto\n") ;
    fprintf(fpo,"      xorig zoomymax tracecad{epcadre add}if lineto\n") ;
    fprintf(fpo,"      closepath} def cadrepath clip\n") ;
    fprintf(fpo,"    }{ nbrcad 4 mul copy\n") ;
    fprintf(fpo,"      /xxmin{dup}def /yymin{dup}def /xxmax{dup}def /yymax{dup}def\n") ;
    fprintf(fpo,"      /t@stx@{ dup yymax ge{/yymax exch def}{pop}ifelse dup yymin le{/yymin exch def}{pop}ifelse\n") ;
    fprintf(fpo,"               dup xxmax ge{/xxmax exch def}{pop}ifelse dup xxmin le{/xxmin exch def}{pop}ifelse}def\n") ;
    fprintf(fpo,"      nbrcad { t@stx@ } repeat 14 taillefonte\n") ;
    fprintf(fpo,"      /xxmin xxmin tfont 1.2 mul sub def /xxmax xxmax tfont 1.2 mul add def\n") ;
    fprintf(fpo,"      /yymin yymin tfont 1.2 mul sub def /yymax yymax tfont 1.2 mul add def\n") ;
    fprintf(fpo,"      yymin zoomymin sub dup vsize ge {/yorig zoomymin vsize add def /xorig zoomxmin def pop}\n") ;
    fprintf(fpo,"      {zoomymax yymax sub dup vsize ge {/yorig zoomymax def /xorig zoomxmin def pop pop}\n") ;
    fprintf(fpo,"      {xxmin zoomxmin sub dup vsize ge {/angle 90 def /xorig zoomxmin def /yorig zoomymin def pop pop pop}\n") ;
    fprintf(fpo,"      {zoomxmax xxmax sub dup vsize ge {/angle 90 def /xorig zoomxmax vsize sub def\n") ;
    fprintf(fpo,"                                        /yorig zoomymax hsize sub def pop pop pop pop}\n") ;
    fprintf(fpo,"      {/angle 90 def /xorig zoomxmax vsize sub def /yorig zoomymax hsize sub def\n") ;
    fprintf(fpo,"      2 copy ge {pop /xorig zoomxmin def /yorig zoomymin def}{exch pop} ifelse\n") ;
    fprintf(fpo,"      2 copy ge {pop /angle 0 def /yorig zoomymax def /xorig zoomxmin def}{exch pop} ifelse\n") ;
    fprintf(fpo,"      ge {/angle 0 def /yorig zoomymin vsize add def /xorig zoomxmin def}if}ifelse}ifelse}ifelse}ifelse\n") ;
    fprintf(fpo,"      gsave xorig yorig moveto angle rotate hsize 0 rlineto 0 vsize neg rlineto hsize neg 0 rlineto closepath\n") ;
    fprintf(fpo,"      1.0 setgray fill grestore\n") ;
    fprintf(fpo,"    } ifelse\n") ;
    fprintf(fpo,"    /sh@w {gsave tfont -0.1 mul dup rmoveto dup stringwidth pop tfont 0.2 mul add dup 0 rlineto\n") ;
    fprintf(fpo,"    0 tfont 1.2 mul rlineto neg 0 rlineto closepath 1.0 setgray fill grestore show}def\n") ;
    fprintf(fpo,"    nbrcad -1 1{5 1 roll\n") ;
    fprintf(fpo,"    /ymax exch 10 mul round 10 div def /ymin exch 10 mul round 10 div def\n") ;
    fprintf(fpo,"    /xmax exch 10 mul round 10 div def /xmin exch 10 mul round 10 div def\n") ;
    fprintf(fpo,"    newpath\n") ;
    fprintf(fpo,"    xmin c@rli 2 div sub ymin c@rli 2 div sub moveto xmax c@rli 2 div add ymin c@rli 2 div sub lineto\n") ;
    fprintf(fpo,"    xmax c@rli 2 div add ymax c@rli 2 div add lineto xmin c@rli 2 div sub ymax c@rli 2 div add lineto\n") ;
    fprintf(fpo,"    closepath stroke\n") ;
    fprintf(fpo,"    20 taillefonte xmin ymax moveto tfont dup 0.1 mul exch -1.1 mul rmoveto\n") ;
    fprintf(fpo,"    /st 2 string def dup st cvs sh@w 14 taillefonte\n") ;
    fprintf(fpo,"    gsave xmin c@rli sub ymin c@rli sub moveto\n") ;
    fprintf(fpo,"    tfont dup 0.1 mul exch -1.1 mul rmoveto /st 10 string def xmin st cvs sh@w grestore\n") ;
    fprintf(fpo,"    gsave xmin c@rli sub ymin c@rli sub translate 90 rotate\n") ;
    fprintf(fpo,"    tfont dup 0.1 mul exch 0.1 mul moveto /st 10 string def ymin st cvs sh@w grestore\n") ;
    fprintf(fpo,"    gsave xmax c@rli add ymax c@rli add moveto\n") ;
    fprintf(fpo,"    /st 10 string def xmax st cvs dup stringwidth pop neg tfont 0.1 mul sub tfont 0.1 mul rmoveto\n") ;
    fprintf(fpo,"    sh@w grestore\n") ;
    fprintf(fpo,"    gsave xmax c@rli add ymax c@rli add translate 90 rotate\n") ;
    fprintf(fpo,"    /st 10 string def ymax st cvs dup stringwidth pop neg tfont 0.1 mul sub tfont -1.1 mul moveto\n") ;
    fprintf(fpo,"    sh@w grestore 17 taillefonte\n") ;
    fprintf(fpo,"    gsave xorig yorig moveto angle rotate dup 1 add tfont mul -1.2 mul 0 exch rmoveto\n") ;
    fprintf(fpo,"    st cvs show ( : ) show xmin st cvs tab xmax st cvs tab ymin st cvs tab ymax st cvs tab grestore\n") ;
    fprintf(fpo,"    } for\n") ;
    fprintf(fpo,"    gsave xorig yorig moveto angle rotate 0 tfont -1.2 mul rmoveto\n") ;
    fprintf(fpo,"    (1 : ) stringwidth rmoveto (xmin) tab (xmax) tab (ymin) tab (ymax) tab\n") ;
    fprintf(fpo,"  } def\n") ;
    fprintf(fpo,"}{ %% macros strandard\n") ;
    fprintf(fpo,"  /m@veto {moveto}  bind def /m@v@to  {moveto}  bind def /m@v@t@ {moveto}  bind def\n") ;
    fprintf(fpo,"  /sh@w   {show}    bind def /l@neto  {lineto}  bind def /l@n@to {lineto}  bind def\n") ;
    fprintf(fpo,"  /rl@neto{rlineto} bind def /rl@n@to {rlineto} bind def /t@stel { }       def\n") ;
    fprintf(fpo,"  /t@st@l { }       def      /autoreca{ }       def      /encadre{ }       def\n") ;
    fprintf(fpo,"} ifelse\n") ;
    fprintf(fpo,"%%%%\n") ;
    fprintf(fpo,"%%%% Procedures de trace\n") ;
    fprintf(fpo,"lartrait 0 eq{ /str@k { } def }{ /str@k { stroke } bind def } ifelse\n") ;
    fprintf(fpo,"couleur not { /fondcou false def } if \n") ;
    fprintf(fpo,"fondcou {\n") ;
    fprintf(fpo,"  /s@tc@u { setrgbcolor } bind def\n") ;
    if ( modepsc != 0 ){
      fprintf(fpo,"  /s@tcou { setrgbcolor } bind def\n") ;
    }
    fprintf(fpo,"  /setc@u { setrgbcolor } bind def\n") ;
    fprintf(fpo,"}{\n") ;
    fprintf(fpo,"  /s@tc@u { %% inversion du blanc en noir et blanc en noir\n") ;
    fprintf(fpo,"    3 copy add add dup 0 eq { pop pop pop pop 1.0 1.0 1.0 }\n") ;
    fprintf(fpo,"                            { 3 eq { pop pop pop 0.0 0.0 0.0 } if } ifelse\n") ;
    fprintf(fpo,"    setrgbcolor\n") ;
    fprintf(fpo,"  } bind def\n") ;
    if ( modepsc != 0 ){
      fprintf(fpo,"  /s@tcou { %% inversion du noir en blanc\n") ;
      fprintf(fpo,"    3 copy add add 0 eq { pop pop pop 1.0 1.0 1.0 } if\n") ;
      fprintf(fpo,"    setrgbcolor\n") ;
      fprintf(fpo,"  } bind def\n") ;
    }
    fprintf(fpo,"  /setc@u { %% inversion du blanc en noir (pour T)\n") ;
    fprintf(fpo,"    3 copy add add 3 eq { pop pop pop 0.0 0.0 0.0 } if\n") ;
    fprintf(fpo,"    setrgbcolor\n") ;
    fprintf(fpo,"  } bind def\n") ;
    fprintf(fpo,"} ifelse\n") ;
    fprintf(fpo,"couleur { %% couleur\n") ;
    fprintf(fpo,"/setcose{ pop s@tc@u } bind def\n") ;
    fprintf(fpo,"/setcofa{ pop s@tc@u } bind def\n") ;
    fprintf(fpo,"}{ %% Noir&Blanc : palette de gris\n") ;
    fprintf(fpo,"/setcose{\n") ;
    fprintf(fpo,"  dup 1 exch sub gsefo mul exch gsecl mul add setgray pop pop pop\n") ;
    fprintf(fpo,"} bind def\n") ;
    fprintf(fpo,"/setcofa{\n") ;
    fprintf(fpo,"  dup 1 exch sub gfafo mul exch gfacl mul add setgray pop pop pop\n") ;
    fprintf(fpo,"} bind def\n") ;
    fprintf(fpo,"} ifelse\n") ;
    if ( iep != 0 ){
      fprintf(fpo,"/epais{ %% epaisseur du trait         pile : num\n") ;
      fprintf(fpo,"  reduc mul setlinewidth\n") ;
      fprintf(fpo,"} bind def\n") ;
    }
    fprintf(fpo,"/S{ %% trace d'une suite de segments couleurs pile : xn yn .. x0 y0 n coul\n") ;
    fprintf(fpo,"  gsave setcose newpath 3 1 roll m@veto { l@neto } repeat str@k grestore\n") ;
    fprintf(fpo,"} bind def\n") ;
    if ( iPo != 0 ){
      fprintf(fpo,"/P{ %% trace d'un polygone ferme en couleur  pile : xn yn .. x0 y0 n coul\n") ;
      fprintf(fpo,"  gsave setcose newpath 3 1 roll m@veto 1 sub { l@neto } repeat closepath str@k grestore\n") ;
      fprintf(fpo,"} bind def\n") ;
    }
    if ( iFa != 0 ){
      fprintf(fpo,"/F{ %% remplissage d'une facette  pile : xn yn .. x1 y1 n coul\n") ;
      fprintf(fpo,"  gsave setcofa newpath 3 1 roll m@v@to 1 sub { l@n@to } repeat closepath\n") ;
      fprintf(fpo,"  tracecou dup 0 ne { 1 eq { 1 setgray } if fill }{ pop } ifelse grestore\n") ;
      fprintf(fpo,"} bind def\n") ;
    }
    if ( iFP != 0 ){
      fprintf(fpo,"/FP{ %% facette et polygone pile : xn yn .. x1 y1 n coul_pol coul_fac\n") ;
      fprintf(fpo,"  gsave setcofa  5 -1 roll dup 6 1 roll 2 mul 5 add 4 roll\n") ;
      fprintf(fpo,"  newpath 3 1 roll m@veto 1 sub { l@neto } repeat closepath\n") ;
      fprintf(fpo,"  gsave tracecou dup 0 ne { 1 eq { 1 setgray } if fill }{ pop } ifelse grestore\n") ;
      fprintf(fpo,"  setcose str@k grestore\n") ;
      fprintf(fpo,"} bind def\n") ;
    }
    if ( ire != 0 ){
      fprintf(fpo,"/r{ %% trace du bord d'un rectangle en coul pile : dx dy x y coul\n") ;
      fprintf(fpo,"  gsave setcose 0 setlinejoin newpath\n") ;
      fprintf(fpo,"  m@veto 0 exch 0 exch 2 copy neg 6 2 roll rl@neto rl@neto rl@neto closepath\n") ;
      fprintf(fpo,"  str@k grestore\n") ;
      fprintf(fpo,"} bind def\n") ;
    }
    if ( iRe != 0 ){
      fprintf(fpo,"/R{ %% remplissage d'un rectangle en coul  pile : dx dy x y coul\n") ;
      fprintf(fpo,"  gsave setcofa newpath\n") ;
      fprintf(fpo,"  m@v@to 0 exch 0 exch 2 copy neg 6 2 roll rl@n@to rl@n@to rl@n@to closepath\n") ;
      fprintf(fpo,"  fill grestore\n") ;
      fprintf(fpo,"} bind def\n") ;
    }
    if ( iel != 0 ){
      fprintf(fpo,"/el{ %% trace d'une ellipse en coul  pile : angle_deb angle_fin largeur hauteur x y coul\n") ;
      fprintf(fpo,"  gsave setcose t@stel /savematrix matrix currentmatrix def\n") ;
      fprintf(fpo,"  newpath translate scale 0 0 1 5 3 roll arc\n") ;
      fprintf(fpo,"  savematrix setmatrix str@k grestore\n") ;
      fprintf(fpo,"} bind def\n") ;
    }
    if ( iEl != 0 ){
      fprintf(fpo,"/El{ %% remplissage d'une ellipse en coul  pile : angle_deb angle_fin largeur hauteur x y coul\n") ;
      fprintf(fpo,"  gsave setcofa t@st@l /savematrix matrix currentmatrix def\n") ;
      fprintf(fpo,"  newpath translate scale 0 0 1 5 3 roll arc\n") ;
      fprintf(fpo,"  savematrix setmatrix fill grestore\n") ;
      fprintf(fpo,"} bind def\n") ;
    }
    if ( ity != 0 ){
      fprintf(fpo,"/typet{ %% type de trait (continu tirete 1 ou 2\n") ;
      fprintf(fpo,"  [[] [4 4] [4 4]] exch get 0 setdash\n") ;
      fprintf(fpo,"} bind def\n") ;
    }
    if ( iTe != 0 ){
      fprintf(fpo,"traceecr { %% trace effectif des chaines de caracteres\n") ;
      fprintf(fpo,"couleur { %% couleur\n") ;
      fprintf(fpo,"  /T{ pop setc@u m@v@t@ sh@w } bind def\n") ;
      fprintf(fpo,"}{\n") ;
      fprintf(fpo,"  /T{ setgray pop pop pop m@v@t@ sh@w } bind def\n") ;
      fprintf(fpo,"} ifelse\n") ;
      fprintf(fpo,"}{ %% Pas de trace des chaines de caracteres\n") ;
      fprintf(fpo,"  /T{ 7 { pop } repeat } bind def\n") ;
      fprintf(fpo,"} ifelse\n") ;
    }
    fprintf(fpo,"/cadrepath{ %% definition du cadre de trace  pile : xmin ymin xmax ymax\n") ;
    fprintf(fpo,"  newpath 4 copy moveto 3 1 roll exch lineto 4 2 roll lineto lineto closepath\n") ;
    fprintf(fpo,"} bind def\n\n") ;
    fprintf(fpo,"%%%%EndProcSet \n\n") ;
    fprintf(fpo,"texval not {  %% orientation du dessin sur A4\n") ;
    fprintf(fpo,"  21.0 cm 2 div 29.7 cm 2 div translate\n") ;
    fprintf(fpo,"  vsize hsize gt {\n") ;
    fprintf(fpo,"    zoomxmax zoomxmin sub zoomymax zoomymin sub gt {\n") ;
    fprintf(fpo,"      90 rotate vsize hsize /vsize exch def /hsize exch def } if\n") ;
    fprintf(fpo,"  }{\n") ;
    fprintf(fpo,"    zoomxmax zoomxmin sub zoomymax zoomymin sub lt {\n") ;
    fprintf(fpo,"      90 rotate vsize hsize /vsize exch def /hsize exch def } if\n") ;
    fprintf(fpo,"  } ifelse\n") ;
    fprintf(fpo,"} if \n") ;
    fprintf(fpo,"/reduc %% calcul du coeff de reduction du a un eventuel zoom\n") ;
    fprintf(fpo,"  vsize zoomymax zoomymin sub div hsize zoomxmax zoomxmin sub div ge {\n") ;
    fprintf(fpo,"    zoomxmax zoomxmin sub %5i div\n",xpixels) ;
    fprintf(fpo,"    vsize %5i div hsize %5i div 2 copy le { div mul }{pop pop} ifelse\n",ypixels,xpixels) ;
    fprintf(fpo,"  }{\n") ;
    fprintf(fpo,"    zoomymax zoomymin sub %5i div\n",ypixels) ;
    fprintf(fpo,"    hsize %5i div vsize %5i div 2 copy le { div mul }{pop pop} ifelse\n",xpixels,ypixels) ;
    fprintf(fpo,"  } ifelse\n") ;
    fprintf(fpo,"def\n") ;
    fprintf(fpo,"lartrait -1 eq {\n") ;
    fprintf(fpo,"  /reduclog reduc def\n") ;
    fprintf(fpo,"  texval { /reduc reduc reductex mul def } if\n") ;
    fprintf(fpo,"}{\n") ;
    fprintf(fpo," /reduclog reduc texval {reductex mul} if def /reduc lartrait def\n") ;
    fprintf(fpo,"} ifelse\n") ;
    fprintf(fpo,"/epcadre epcadre reduc mul def\n") ;
    fprintf(fpo,"texval not\n") ;
    fprintf(fpo,"{ %% sortie directe sur imprimante\n") ;
    fprintf(fpo,"  vsize zoomymax zoomymin sub epcadre 2 mul add div hsize zoomxmax zoomxmin sub epcadre 2 mul add div\n") ;
    fprintf(fpo,"  2 copy gt { exch } if pop dup scale\n") ;
    fprintf(fpo,"  zoomxmin zoomxmax add 2 div neg zoomymin zoomymax add 2 div neg translate\n") ;
    fprintf(fpo,"}\n") ;
    fprintf(fpo,"{ %% integration sur TeX\n") ;
    fprintf(fpo,"  /bavard false def\n") ;      
    fprintf(fpo,"  texdvipsval\n") ;
    fprintf(fpo,"  { %% integration par DVIPS et DVI2PS\n") ;
    fprintf(fpo,"    -1 -1 scale\n") ;
    fprintf(fpo,"    zoomxmax neg tracecad {epcadre sub} if zoomymin neg tracecad {epcadre add} if translate\n") ;
    fprintf(fpo,"  }\n") ;
    fprintf(fpo,"  { %% integration par TEXTURES\n") ;
    fprintf(fpo,"    zoomxmin neg tracecad {epcadre add} if zoomymin neg tracecad {epcadre add} if translate\n") ;
    fprintf(fpo,"  }\n") ;
    fprintf(fpo,"  ifelse\n") ;
    fprintf(fpo,"}\n") ;
    fprintf(fpo,"ifelse\n") ;
    fprintf(fpo,"gsave\n") ;
    fprintf(fpo,"zoomxmin zoomymin zoomxmax zoomymax cadrepath clip newpath\n") ;
    fprintf(fpo,"/charge { %% Chargement d'une fonte. pile : /police /corps (/ si normal) /style (/ si roman) haut larg (.) (mono ou prop)\n") ;
    fprintf(fpo,"  /mat [10 0 0 10 0 0] def /mats [10 0 0 10 0 0] def 6 3 roll\n") ;
    fprintf(fpo,"  /g@t {getinterval} bind def /p@t {putinterval} bind def\n") ;
    fprintf(fpo,"  /sty exch 20 string cvs def /cor exch 20 string cvs def /pol exch 30 string cvs def\n") ;
    fprintf(fpo,"  /lpo pol length def /lco cor length def /lst sty length def /fon 71 string def\n") ;
    fprintf(fpo,"  fon 0 pol p@t lco lst add 0 ne\n") ;
    fprintf(fpo,"   {fon lpo (-) p@t fon lpo 1 add cor p@t fon lpo 1 add lco add sty p@t}{/lpo lpo 1 sub def} ifelse\n") ;
    fprintf(fpo,"  bavard{/err 100 string def err 0 fon 0 lpo lco lst 1 add add add g@t p@t} if\n") ;
    fprintf(fpo,"  lst 0 ne{sty 0 1 g@t (O) eq {mats 2 mats 3 get 0.25 mul put} if\n") ;
    fprintf(fpo,"           sty 0 1 g@t (I) eq {mats 2 mats 3 get 0.2  mul put} if} if\n") ;
    fprintf(fpo,"  fon 0 lpo lco lst 1 add add add g@t cvn FontDirectory exch known\n") ;
    fprintf(fpo,"    {/lfo lpo lco lst 1 add add add def bavard {err 50 (OK) p@t err counttomark 1 roll} if}\n") ;
    fprintf(fpo,"    {fon 0 lpo lco 1 add add g@t cvn FontDirectory exch known {true}{\n") ;
    fprintf(fpo,"      fon 0 lpo g@t cvn FontDirectory exch known {/cor () def /lco 0 def true}\n") ;
    fprintf(fpo,"      {/cor (Roman) def /lco 5 def fon 0 fon lpo 1 add cor p@t\n") ;
    fprintf(fpo,"       lpo lco 1 add add g@t cvn FontDirectory exch known\n") ;
    fprintf(fpo,"      {true}{/cor () def /lco 0 def false} ifelse } ifelse } ifelse\n") ;
    fprintf(fpo,"      {fon lpo (-) p@t fon lpo 1 add cor p@t fon lpo 1 add lco add sty p@t\n") ;
    fprintf(fpo,"        lst 0 ne {fon 0 lpo lco lst 1 add add add g@t cvn FontDirectory exch known}{false}ifelse\n") ;
    fprintf(fpo,"        {/lfo lpo lco lst 1 add add add def\n") ;
    fprintf(fpo,"          bavard {err 50 (--> ) p@t err 54 fon 0 lfo g@t p@t err counttomark 1 roll } if }\n") ;
    fprintf(fpo,"        {lco 0 ne {/lfo lpo lco 1 add add def}{/lfo lpo def} ifelse /mat mats def\n") ;
    fprintf(fpo,"          bavard {err 50 (--> ) p@t err 54 fon 0 lfo g@t p@t\n") ;
    fprintf(fpo,"                  lst 0 ne{err 54 lfo add ( /penche/) p@t} if err counttomark 1 roll} if } ifelse }\n") ;
    fprintf(fpo,"      {dup dup /fon exch (m) eq {fontm0}{fontp0}ifelse 71 string cvs def /lfo fon length def /mat mats def\n") ;
    fprintf(fpo,"         bavard {(m) eq {err 50 (--> fontm0 par defaut) p@t}{err 50 (--> fontp0 par defaut) p@t}ifelse\n") ;
    fprintf(fpo,"                lst 0 ne{err 71 ( /penche/) p@t} if err counttomark 1 roll}{pop} ifelse } ifelse\n") ;
    fprintf(fpo,"    } ifelse\n") ;
    fprintf(fpo,"    pop fon 0 lfo g@t cvn findfont mat makefont setfont newpath 0 0 moveto\n") ;
    fprintf(fpo,"    (abcdefghijklmnopqrstuvwxyz ABCDEFGHIJKLMNOPQRSTUVWXYZ 1234567890) true charpath flattenpath pathbbox\n") ;
    fprintf(fpo,"    4 1 roll exch pop exch sub exch 4 1 roll exch 1.0 mul exch div 3 1 roll 1.0 mul dup /tfont exch def exch div\n") ;
    fprintf(fpo,"    dup mat 2 mat 2 get 4 -1 roll mul put mat 3 mat 3 get 4 -1 roll mul put mat 0 mat 0 get 4 -1 roll mul put\n") ;
    fprintf(fpo,"    fon 0 lfo g@t cvn findfont mat makefont setfont\n") ;
    fprintf(fpo,"} bind def\n") ;
    fprintf(fpo,"%s",fontcour) ;
    fprintf(fpo,"fondcou { gsave zoomxmin zoomymin zoomxmax zoomymax cadrepath 0 setgray fill grestore } if\n") ;
    fprintf(fpo,"0 setlinecap\n") ;
    fprintf(fpo,"1 setlinejoin\n") ;
    do { ligne = fgets(buf,255,fpi) ;
      if ( ligne == buf ){
        fprintf(fpo,"%s",buf) ;
      }
    } while ( ligne == buf ) ;
    fclose(fpi) ;
    fpi = NULL ;

    if ( modepsc != 0 ){
    /* ecriture des legendes qualites ou isovaleurs */
    if ( strlen(chaine[7]) != 0 ){
      fprintf(fpo,"\n%%%% Debut : Legende des qualites ou des isovaleurs\n\n") ;
      fprintf(fpo,"tracemen 1 eq tracecou 2 eq and { %% test sur le trace effectif\n") ;
      fprintf(fpo,"%%%% Redefinition locale pour la legende des qualites ou des isovaleurs\n") ;
      fprintf(fpo,"/epais{ %% epaisseur du trait         pile : num\n") ;
      fprintf(fpo,"  setlinewidth\n") ;
      fprintf(fpo,"} bind def\n") ;
      fprintf(fpo,"/S{ %% trace d'un segment         pile : x2 y2 x1 y1 coul\n") ;
      fprintf(fpo,"  gsave setcose newpath m@veto l@neto stroke grestore\n") ;
      fprintf(fpo,"} bind def\n") ;
      fprintf(fpo,"/R{ %% remplissage d'un rectangle en coul  pile : dx dy x y coul\n") ;
      fprintf(fpo,"  gsave couleur { pop s@tcou\n") ;
      fprintf(fpo,"  }{ dup 1 exch sub gfafo mul exch gfacl mul add setgray pop pop pop\n") ;
      fprintf(fpo,"  } ifelse\n") ;
      fprintf(fpo,"  newpath m@veto 0 exch 0 exch 2 copy neg 6 2 roll rl@neto rl@neto rl@neto\n") ;
      fprintf(fpo,"  closepath gsave fill grestore 0 setgray 0.5 epais stroke grestore\n") ;
      fprintf(fpo,"} bind def\n") ;
      fprintf(fpo,"/T{ %% impression chaine de caracteres\n") ;
      fprintf(fpo,"  couleur { pop setrgbcolor }{ setgray pop pop pop } ifelse\n") ;
      fprintf(fpo,"  m@v@t@ sh@w\n") ;
      fprintf(fpo,"} bind def\n") ;
      fprintf(fpo,"0 setlinecap\n") ;
      fprintf(fpo,"0 setlinejoin\n\n") ;
      fprintf(fpo,"%%%% Legende\n") ;
      fprintf(fpo,"%s",chaine[7]) ;          
      fprintf(fpo,"} if\n") ;
      fprintf(fpo,"\n%%%% Fin   : Legende des qualites ou des isovaleurs\n") ;          
    }

    /* ecriture des menus */
    if ( menu ){
      fprintf(fpo,"\n%%%% Debut : Menus\n\n") ;          
      fprintf(fpo,"tracemen 2 eq { %% test sur le trace effectif\n") ;
      fprintf(fpo,"%%%% Redefinition locale pour les menus\n") ;          
      fprintf(fpo,"/epais{ %% epaisseur du trait         pile : num\n") ;
      fprintf(fpo,"  dup 4 ge { 0.2 setgray }{ 0 setgray } ifelse setlinewidth\n") ;
      fprintf(fpo,"} bind def\n") ;
      fprintf(fpo,"/S{ %% trace d'un segment         pile : x2 y2 x1 y1 coul\n") ;
      fprintf(fpo,"  gsave couleur { pop s@tc@u }{ pop pop pop pop } ifelse\n") ;
      fprintf(fpo,"  newpath moveto lineto stroke grestore\n") ;
      fprintf(fpo,"} bind def\n") ;
      fprintf(fpo,"/R{ %% remplissage d'un rectangle en coul  pile : dx dy x y coul\n") ;
      fprintf(fpo,"  gsave couleur { pop setrgbcolor }\n") ;
      fprintf(fpo,"  { dup 1 exch sub gfafo mul exch gfacl mul add setgray pop pop pop\n") ;
      fprintf(fpo,"  } ifelse\n") ;
      fprintf(fpo,"  newpath moveto 0 exch 0 exch 2 copy neg 6 2 roll rl@neto rl@neto rl@neto\n") ;
      fprintf(fpo,"  closepath gsave fill grestore 0 setgray 0.5 epais stroke grestore\n") ;
      fprintf(fpo,"} bind def\n") ;
      fprintf(fpo,"/T{ %% impression chaine de caracteres\n") ;
      fprintf(fpo,"  couleur { pop setrgbcolor }{ setgray pop pop pop } ifelse\n") ;
      fprintf(fpo,"  m@v@t@ sh@w\n") ;
      fprintf(fpo,"} bind def\n") ;
      fprintf(fpo,"0 setlinecap\n") ;
      fprintf(fpo,"0 setlinejoin\n\n") ;
      car = "(" ;
      for ( i = 0 ; i < strlen(chaine[5]) ; i++ ){
        if ( chaine[5][i] == *car ){
          sprintf(&chaine[5][i],"%s%s",
                 "(?                                                                                              .)",
                  &chaine[5][i+98]);
        }             
      }
      for ( i = 0 ; i < 7 ; i++ ){
        if ( strlen(chaine[i]) != 0 ){
          switch (i){
            case 0 :
              sprintf(format,"Documentation") ;
              break ;
            case 1 :
              sprintf(format,"Historique") ;
              break ;
            case 2 :
              sprintf(format,"Menu") ;
              break ;
            case 3 :
              sprintf(format,"Lignes lues") ;
              break ;
            case 4 :
              sprintf(format,"Invite") ;
              break ;
            case 5 :
              sprintf(format,"Ligne de saisie") ;
              break ;
            case 6 :
              sprintf(format,"Erreur") ;
              break ;
          }
          fprintf(fpo,"\n%%%% Debut : %s\n",format) ;          
          fprintf(fpo,"%s",chaine[i]) ;          
          fprintf(fpo,"%%%% Fin   : %s\n",format) ;          
        }  
      }
      fprintf(fpo,"} if\n") ;
    }
    }
    fprintf(fpo,"%%\n") ;
    fprintf(fpo,"grestore\n") ;
    fprintf(fpo,"tracecad { %% then trace d'un cadre\n") ;
    fprintf(fpo,"  gsave zoomxmin epcadre sub zoomymin epcadre sub\n") ;
    fprintf(fpo,"  zoomxmax epcadre add zoomymax epcadre add cadrepath clip\n") ;
    fprintf(fpo,"  zoomxmin epcadre 2 div sub zoomymin epcadre 2 div sub\n") ;
    fprintf(fpo,"  zoomxmax epcadre 2 div add zoomymax epcadre 2 div add cadrepath\n") ;
    fprintf(fpo,"  epcadre setlinewidth str@k grestore\n") ;
    fprintf(fpo,"} if\n") ;
    if ( modepsc != 0 ){
      fprintf(fpo,"tracelog { %% then trace du logo\n") ;
      fprintf(fpo,"  gsave /taille 45 reduclog mul def /ZapfChancery-MediumItalic findfont taille scalefont setfont\n") ;
      fprintf(fpo,"  /tailbl taille 16 div def /larlog (Me\\256sto) stringwidth pop tailbl 6 mul add def\n") ;
      fprintf(fpo,"  /haulog taille tailbl 3 mul add def\n") ;
      fprintf(fpo,"  /rays{ gsave larlog 2 div haulog neg translate\n") ;
      fprintf(fpo,"    newpath 121 {larlog 2 mul 0 moveto 0 0 lineto 0.5 rotate larlog 2 mul 0 lineto closepath\n") ;
      fprintf(fpo,"    fill 0 0  moveto 1.0 rotate } repeat grestore\n") ;
      fprintf(fpo,"  } bind def\n") ;
      fprintf(fpo,"  zoomxmax zoomxmin add larlog sub 2 div zoomxmax zoomxmin sub larlog sub 2 div xpolog mul add\n") ;
      fprintf(fpo,"  zoomymax zoomymin add haulog add 2 div zoomymax zoomymin sub haulog sub 2 div ypolog mul add\n") ;
      fprintf(fpo,"  translate 0 0 moveto larlog haulog neg 0 exch 0 exch 2 copy neg 6 2 roll\n") ;
      fprintf(fpo,"  rlineto rlineto rlineto closepath gsave fondcou {0}{1} ifelse setgray fill grestore\n") ;
      fprintf(fpo,"  1 reduclog mul setlinewidth fondcou {1}{0} ifelse setgray stroke\n") ;
      fprintf(fpo,"  tailbl 2 mul taille 0.75 mul tailbl add neg translate 0 0 moveto\n") ;
      fprintf(fpo,"  gsave (Me\\256sto) true charpath clip rays grestore gsave (Me) stringwidth pop 0 rmoveto\n") ;
      fprintf(fpo,"  (e) stringwidth pop (\\302) stringwidth pop add 2 div neg 0 rmoveto\n") ;
      fprintf(fpo,"  (\\302) true charpath clip rays grestore grestore\n") ;
      fprintf(fpo,"} if\n") ;
    }
    fprintf(fpo,"encadre\n") ;
    fprintf(fpo,"gsave grestore\n") ;
    fprintf(fpo,"texval not { showpage } if\n") ;
    fprintf(fpo,"bavard{\n") ;
    fprintf(fpo,"  1 cm 28.7 cm translate /Courier findfont 10 scalefont setfont\n") ;
    fprintf(fpo,"  0 -15 moveto (FONTE DEMANDEE                                        FONTE UTILISEE) show\n") ;
    fprintf(fpo,"  0 -30 moveto counttomark {dup show stringwidth pop neg -15 rmoveto} repeat\n") ;
    fprintf(fpo,"  showpage\n") ;
    fprintf(fpo,"} if\n") ;
    fprintf(fpo,"cleartomark\n") ;
    fclose(fpo) ;
    fpo = NULL ;
    if ( modepsc != 0 ){
      for (i = 0; i < 8 ; i++) { chaine[i] = '\0' ; free(chaine[i]) ; }
    }
    remove("TEMPORAIRE.EPS") ;
    menu = 0 ;
    /* sauvegarde du fichier temporaire des histogrammes de qualite */
    if ( modepsc != 0 ){
      if ((fpi = fopen("TEMPORAIRE.QUA","r"))!=NULL) {
        sprintf(format,"%%.%ds.qua",*length) ;
        sprintf(buf,format,nomfichier) ;
        fpo = NULL ;
        if ((fpo = fopen(buf,"w"))==NULL) {
          printf ("Erreur d'ouverture de fichier\n") ; 
          *length = -1 ;
        } else {
          do { ligne = fgets(buf,255,fpi) ;
            if ( ligne == buf ){
              fprintf(fpo,"%s",buf) ;
            }
          } while ( ligne == buf ) ;
          fclose(fpi) ;
          fpi = NULL ;
          fclose(fpo) ;
          fpo = NULL ; 
          *length = 0 ;
          remove("TEMPORAIRE.QUA") ;
        }
      } else {
        *length = 0 ;
      }
    }
  }
}
