/*    definition de la macro call
/*    elle sert quand des procedures c et fortran sont melangees
*/

/* definition vieillotte de la concatenation
/* dans l'ANSI C on peut faite mieux */
#ifdef __STDC__
#define name2(a,b)  a##b
#else
#ifdef _SILICO_
#define name2(a,b)  a ## b
#else
#ifdef BSD  /*BSD way: ok pour dn10000*/
#define name2(a,b) a\
b

#else /*System V way:*/
#define name2(a,b) a/**/b
#endif
#endif
#endif


/*  SI VOTRE COMPILATEUR FORTRAN PREFERE GENERE UN _ APRES */
/*   LES NOMS DES SUBROUTINES, ALORS PRENEZ CETTE DEFINITION DE call */
/*  #define call(x)  name2(x,_)   */
/*  SINON PRENEZ CELLE CI  */
/* #define call(x)  x  */

#ifdef F77_NO_UNDER_SCORE
#define call(x) x
#else
#define call(x) name2(x,_)
#endif
