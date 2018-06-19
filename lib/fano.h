#define	LL 1	/* Select Layland-Lushbaugh code */

/* Convolutional coding polynomials. All are rate 1/2, K=32 */
#ifdef	NASA_STANDARD
/* "NASA standard" code by Massey & Costello
 * Nonsystematic, quick look-in, dmin=11, dfree=23
 * used on Pioneer 10-12, Helios A,B
 */
#define	POLY1	0xbbef6bb7
#define	POLY2	0xbbef6bb5
#endif

#ifdef	MJ
/* Massey-Johannesson code
 * Nonsystematic, quick look-in, dmin=13, dfree>=23
 * Purported to be more computationally efficient than Massey-Costello
 */
#define	POLY1	0xb840a20f
#define POLY2	0xb840a20d
#endif

#ifdef	LL
/* Layland-Lushbaugh code
 * Nonsystematic, non-quick look-in, dmin=?, dfree=?
 */
#define	POLY1	0xf2d05351
#define	POLY2	0xe4613c47
#endif


