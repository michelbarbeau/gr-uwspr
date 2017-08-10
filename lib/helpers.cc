/* -*- c++ -*- */
/*
 * Copyright 2017 Michel Barbeau, Carleton University.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#include "helpers.h"
#include <stdio.h>
#include <ctype.h>

namespace gr {
   namespace uwspr {

    #define HASH_LITTLE_ENDIAN 1

    #define hashsize(n) ((uint32_t)1<<(n))
    #define hashmask(n) (hashsize(n)-1)
    #define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))

    /*
       -------------------------------------------------------------------------------
       mix -- mix 3 32-bit values reversibly.

       This is reversible, so any information in (a,b,c) before mix() is
       still in (a,b,c) after mix().

       If four pairs of (a,b,c) inputs are run through mix(), or through
       mix() in reverse, there are at least 32 bits of the output that
       are sometimes the same for one pair and different for another pair.
       This was tested for:
     * pairs that differed by one bit, by two bits, in any combination
       of top bits of (a,b,c), or in any combination of bottom bits of
       (a,b,c).
     * "differ" is defined as +, -, ^, or ~^.  For + and -, I transformed
       the output delta to a Gray code (a^(a>>1)) so a string of 1's (as
       is commonly produced by subtraction) look like a single 1-bit
       difference.
     * the base values were pseudorandom, all zero but one bit set, or
       all zero plus a counter that starts at zero.

       Some k values for my "a-=c; a^=rot(c,k); c+=b;" arrangement that
       satisfy this are
        4  6  8 16 19  4
        9 15  3 18 27 15
       14  9  3  7 17  3
       Well, "9 15 3 18 27 15" didn't quite get 32 bits diffing
       for "differ" defined as + with a one-bit base and a two-bit delta.  I
       used http://burtleburtle.net/bob/hash/avalanche.html to choose
       the operations, constants, and arrangements of the variables.

       This does not achieve avalanche.  There are input bits of (a,b,c)
       that fail to affect some output bits of (a,b,c), especially of a.  The
       most thoroughly mixed value is c, but it doesn't really even achieve
       avalanche in c.

       This allows some parallelism.  Read-after-writes are good at doubling
       the number of bits affected, so the goal of mixing pulls in the opposite
       direction as the goal of parallelism.  I did what I could.  Rotates
       seem to cost as much as shifts on every machine I could lay my hands
       on, and rotates are much kinder to the top and bottom bits, so I used
       rotates.
       -------------------------------------------------------------------------------
     */
    #define mix(a,b,c) \
       { \
          a -= c;  a ^= rot(c, 4);  c += b; \
          b -= a;  b ^= rot(a, 6);  a += c; \
          c -= b;  c ^= rot(b, 8);  b += a; \
          a -= c;  a ^= rot(c,16);  c += b; \
          b -= a;  b ^= rot(a,19);  a += c; \
          c -= b;  c ^= rot(b, 4);  b += a; \
       }

    /*
       -------------------------------------------------------------------------------
       final -- final mixing of 3 32-bit values (a,b,c) into c

       Pairs of (a,b,c) values differing in only a few bits will usually
       produce values of c that look totally different.  This was tested for
     * pairs that differed by one bit, by two bits, in any combination
       of top bits of (a,b,c), or in any combination of bottom bits of
       (a,b,c).
     * "differ" is defined as +, -, ^, or ~^.  For + and -, I transformed
       the output delta to a Gray code (a^(a>>1)) so a string of 1's (as
       is commonly produced by subtraction) look like a single 1-bit
       difference.
     * the base values were pseudorandom, all zero but one bit set, or
       all zero plus a counter that starts at zero.

       These constants passed:
       14 11 25 16 4 14 24
       12 14 25 16 4 14 24
       and these came close:
       4  8 15 26 3 22 24
       10  8 15 26 3 22 24
       11  8 15 26 3 22 24
       -------------------------------------------------------------------------------
     */
    #define final(a,b,c) \
       { \
          c ^= b; c -= rot(b,14); \
          a ^= c; a -= rot(c,11); \
          b ^= a; b -= rot(a,25); \
          c ^= b; c -= rot(b,16); \
          a ^= c; a -= rot(c,4);  \
          b ^= a; b -= rot(a,14); \
          c ^= b; c -= rot(b,24); \
       }

    /*
       -------------------------------------------------------------------------------
       hashlittle() -- hash a variable-length key into a 32-bit value
       k       : the key (the unaligned variable-length array of bytes)
       length  : the length of the key, counting by bytes
       initval : can be any 4-byte value
       Returns a 32-bit value.  Every bit of the key affects every bit of
       the return value.  Two keys differing by one or two bits will have
       totally different hash values.

       The best hash table sizes are powers of 2.  There is no need to do
       mod a prime (mod is sooo slow!).  If you need less than 32 bits,
       use a bitmask.  For example, if you need only 10 bits, do
       h = (h & hashmask(10));
       In which case, the hash table should have hashsize(10) elements.

       If you are hashing n strings (uint8_t **)k, do it like this:
       for (i=0, h=0; i<n; ++i) h = hashlittle( k[i], len[i], h);

       By Bob Jenkins, 2006.  bob_jenkins@burtleburtle.net.  You may use this
       code any way you wish, private, educational, or commercial.  It's free.

       Use for hash table lookup, or anything where one collision in 2^^32 is
       acceptable.  Do NOT use for cryptographic purposes.
       -------------------------------------------------------------------------------
     */

    uint32_t helpers::nhash( const void *key, size_t length, uint32_t initval)
    {
       uint32_t a,b,c;                                         /* internal state */
       union { const void *ptr; size_t i; } u;      /* needed for Mac Powerbook G4 */

       /* Set up the internal state */
       a = b = c = 0xdeadbeef + ((uint32_t)length) + initval;

       u.ptr = key;
       if (HASH_LITTLE_ENDIAN && ((u.i & 0x3) == 0)) {
          const uint32_t *k = (const uint32_t *)key;           /* read 32-bit chunks */

          /*------ all but last block: aligned reads and affect 32 bits of (a,b,c) */
          while (length > 12)
          {
             a += k[0];
             b += k[1];
             c += k[2];
             mix(a,b,c);
             length -= 12;
             k += 3;
          }

          /*----------------------------- handle the last (probably partial) block */
          /*
           * "k[2]&0xffffff" actually reads beyond the end of the string, but
           * then masks off the part it's not allowed to read.  Because the
           * string is aligned, the masked-off tail is in the same word as the
           * rest of the string.  Every machine with memory protection I've seen
           * does it on word boundaries, so is OK with this.  But VALGRIND will
           * still catch it and complain.  The masking trick does make the hash
           * noticably faster for short strings (like English words).
           */
    #ifndef VALGRIND

          switch(length)
          {
          case 12: c+=k[2]; b+=k[1]; a+=k[0]; break;
          case 11: c+=k[2]&0xffffff; b+=k[1]; a+=k[0]; break;
          case 10: c+=k[2]&0xffff; b+=k[1]; a+=k[0]; break;
          case 9: c+=k[2]&0xff; b+=k[1]; a+=k[0]; break;
          case 8: b+=k[1]; a+=k[0]; break;
          case 7: b+=k[1]&0xffffff; a+=k[0]; break;
          case 6: b+=k[1]&0xffff; a+=k[0]; break;
          case 5: b+=k[1]&0xff; a+=k[0]; break;
          case 4: a+=k[0]; break;
          case 3: a+=k[0]&0xffffff; break;
          case 2: a+=k[0]&0xffff; break;
          case 1: a+=k[0]&0xff; break;
          case 0: return c;             /* zero length strings require no mixing */
          }

    #else /* make valgrind happy */

          k8 = (const uint8_t *)k;
          switch(length)
          {
          case 12: c+=k[2]; b+=k[1]; a+=k[0]; break;
          case 11: c+=((uint32_t)k8[10])<<16;           /* fall through */
          case 10: c+=((uint32_t)k8[9])<<8;           /* fall through */
          case 9: c+=k8[8];                  /* fall through */
          case 8: b+=k[1]; a+=k[0]; break;
          case 7: b+=((uint32_t)k8[6])<<16;           /* fall through */
          case 6: b+=((uint32_t)k8[5])<<8;           /* fall through */
          case 5: b+=k8[4];                  /* fall through */
          case 4: a+=k[0]; break;
          case 3: a+=((uint32_t)k8[2])<<16;           /* fall through */
          case 2: a+=((uint32_t)k8[1])<<8;           /* fall through */
          case 1: a+=k8[0]; break;
          case 0: return c;
          }

    #endif /* !valgrind */

       } else if (HASH_LITTLE_ENDIAN && ((u.i & 0x1) == 0)) {
          const uint16_t *k = (const uint16_t *)key;           /* read 16-bit chunks */
          const uint8_t  *k8;

          /*--------------- all but last block: aligned reads and different mixing */
          while (length > 12)
          {
             a += k[0] + (((uint32_t)k[1])<<16);
             b += k[2] + (((uint32_t)k[3])<<16);
             c += k[4] + (((uint32_t)k[5])<<16);
             mix(a,b,c);
             length -= 12;
             k += 6;
          }

          /*----------------------------- handle the last (probably partial) block */
          k8 = (const uint8_t *)k;
          switch(length)
          {
          case 12: c+=k[4]+(((uint32_t)k[5])<<16);
             b+=k[2]+(((uint32_t)k[3])<<16);
             a+=k[0]+(((uint32_t)k[1])<<16);
             break;
          case 11: c+=((uint32_t)k8[10])<<16;           /* fall through */
          case 10: c+=k[4];
             b+=k[2]+(((uint32_t)k[3])<<16);
             a+=k[0]+(((uint32_t)k[1])<<16);
             break;
          case 9: c+=k8[8];                     /* fall through */
          case 8: b+=k[2]+(((uint32_t)k[3])<<16);
             a+=k[0]+(((uint32_t)k[1])<<16);
             break;
          case 7: b+=((uint32_t)k8[6])<<16;           /* fall through */
          case 6: b+=k[2];
             a+=k[0]+(((uint32_t)k[1])<<16);
             break;
          case 5: b+=k8[4];                     /* fall through */
          case 4: a+=k[0]+(((uint32_t)k[1])<<16);
             break;
          case 3: a+=((uint32_t)k8[2])<<16;           /* fall through */
          case 2: a+=k[0];
             break;
          case 1: a+=k8[0];
             break;
          case 0: return c;                    /* zero length requires no mixing */
          }

       } else {                       /* need to read the key one byte at a time */
          const uint8_t *k = (const uint8_t *)key;

          /*--------------- all but the last block: affect some 32 bits of (a,b,c) */
          while (length > 12)
          {
             a += k[0];
             a += ((uint32_t)k[1])<<8;
             a += ((uint32_t)k[2])<<16;
             a += ((uint32_t)k[3])<<24;
             b += k[4];
             b += ((uint32_t)k[5])<<8;
             b += ((uint32_t)k[6])<<16;
             b += ((uint32_t)k[7])<<24;
             c += k[8];
             c += ((uint32_t)k[9])<<8;
             c += ((uint32_t)k[10])<<16;
             c += ((uint32_t)k[11])<<24;
             mix(a,b,c);
             length -= 12;
             k += 12;
          }

          /*-------------------------------- last block: affect all 32 bits of (c) */
          switch(length)                 /* all the case statements fall through */
          {
          case 12: c+=((uint32_t)k[11])<<24;
          case 11: c+=((uint32_t)k[10])<<16;
          case 10: c+=((uint32_t)k[9])<<8;
          case 9: c+=k[8];
          case 8: b+=((uint32_t)k[7])<<24;
          case 7: b+=((uint32_t)k[6])<<16;
          case 6: b+=((uint32_t)k[5])<<8;
          case 5: b+=k[4];
          case 4: a+=((uint32_t)k[3])<<24;
          case 3: a+=((uint32_t)k[2])<<16;
          case 2: a+=((uint32_t)k[1])<<8;
          case 1: a+=k[0];
             break;
          case 0: return c;
          }
       }

       final(a,b,c);
       c=(32767&c);

       return c;
    }

    void helpers::unpack50(  char *dat, int32_t *n1, int32_t *n2 )
    {
       int32_t i,i4;

       i=dat[0];
       i4=i&255;
       *n1=i4<<20;

       i=dat[1];
       i4=i&255;
       *n1=*n1+(i4<<12);

       i=dat[2];
       i4=i&255;
       *n1=*n1+(i4<<4);

       i=dat[3];
       i4=i&255;
       *n1=*n1+((i4>>4)&15);
       *n2=(i4&15)<<18;

       i=dat[4];
       i4=i&255;
       *n2=*n2+(i4<<10);

       i=dat[5];
       i4=i&255;
       *n2=*n2+(i4<<2);

       i=dat[6];
       i4=i&255;
       *n2=*n2+((i4>>6)&3);
    }

    int helpers::unpackcall( int32_t ncall, char *call )
    {
       char c[]={'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E',
                 'F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T',
                 'U','V','W','X','Y','Z',' '};
       int32_t n;
       int i;
       char tmp[7];

       n=ncall;
       strcpy(call,"......");
       if (n < 262177560 ) {
          i=n%27+10;
          tmp[5]=c[i];
          n=n/27;
          i=n%27+10;
          tmp[4]=c[i];
          n=n/27;
          i=n%27+10;
          tmp[3]=c[i];
          n=n/27;
          i=n%10;
          tmp[2]=c[i];
          n=n/10;
          i=n%36;
          tmp[1]=c[i];
          n=n/36;
          i=n;
          tmp[0]=c[i];
          tmp[6]='\0';
          // remove leading whitespace
          for(i=0; i<5; i++) {
             if( tmp[i] != c[36] )
                break;
          }
          sprintf(call,"%-6s",&tmp[i]);
          // remove trailing whitespace
          for(i=0; i<6; i++) {
             if( call[i] == c[36] ) {
                call[i]='\0';
             }
          }
       } else {
          return 0;
       }
       return 1;
    }

    int helpers::unpackgrid( int32_t ngrid, char *grid)
    {
       char c[]={'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E',
                 'F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T',
                 'U','V','W','X','Y','Z',' '};
       int dlat, dlong;

       ngrid=ngrid>>7;
       if( ngrid < 32400 ) {
          dlat=(ngrid%180)-90;
          dlong=(ngrid/180)*2 - 180 + 2;
          if( dlong < -180 )
             dlong=dlong+360;
          if( dlong > 180 )
             dlong=dlong+360;
          int nlong = 60.0*(180.0-dlong)/5.0;
          int n1 = nlong/240;
          int n2 = (nlong - 240*n1)/24;
          grid[0] = c[10+n1];
          grid[2]=  c[n2];

          int nlat = 60.0*(dlat+90)/2.5;
          n1 = nlat/240;
          n2 = (nlat-240*n1)/24;
          grid[1]=c[10+n1];
          grid[3]=c[n2];
       } else {
          strcpy(grid,"XXXX");
          return 0;
       }
       return 1;
    }

    int helpers::unpackpfx( int32_t nprefix, char *call)
    {
       char nc, pfx[4]={'\0'}, tmpcall[7];
       int i;
       int32_t n;

       strcpy(tmpcall,call);
       if( nprefix < 60000 ) {
          // add a prefix of 1 to 3 characters
          n=nprefix;
          for (i=2; i>=0; i--) {
             nc=n%37;
             if( (nc >= 0) & (nc <= 9) ) {
                pfx[i]=nc+48;
             }
             else if( (nc >= 10) & (nc <= 35) ) {
                pfx[i]=nc+55;
             }
             else {
                pfx[i]=' ';
             }
             n=n/37;
          }

          char * p = strrchr(pfx,' ');
          strcpy(call, p ? p + 1 : pfx);
          strncat(call,"/",1);
          strncat(call,tmpcall,strlen(tmpcall));

       } else {
          // add a suffix of 1 or 2 characters
          nc=nprefix-60000;
          if( (nc >= 0) & (nc <= 9) ) {
             pfx[0]=nc+48;
             strcpy(call,tmpcall);
             strncat(call,"/",1);
             strncat(call,pfx,1);
          }
          else if( (nc >= 10) & (nc <= 35) ) {
             pfx[0]=nc+55;
             strcpy(call,tmpcall);
             strncat(call,"/",1);
             strncat(call,pfx,1);
          }
          else if( (nc >= 36) & (nc <= 125) ) {
             pfx[0]=(nc-26)/10+48;
             pfx[1]=(nc-26)%10+48;
             strcpy(call,tmpcall);
             strncat(call,"/",1);
             strncat(call,pfx,2);
          }
          else {
             return 0;
          }
       }
       return 1;
    }

    int helpers::unpk_( char *message, char *hashtab, char *call_loc_pow, char *callsign)
    {
       int n1,n2,n3,ndbm,ihash,nadd,noprint=0;
       char grid[5],grid6[7],cdbm[3];

       unpack50(message,&n1,&n2);
       if( !unpackcall(n1,callsign) ) return 1;
       if( !unpackgrid(n2, grid) ) return 1;
       int ntype = (n2&127) - 64;
       callsign[12]=0;
       grid[4]=0;

       /*
          Based on the value of ntype, decide whether this is a Type 1, 2, or
          3 message.

        * Type 1: 6 digit call, grid, power - ntype is positive and is a member
          of the set {0,3,7,10,13,17,20...60}

        * Type 2: extended callsign, power - ntype is positive but not
          a member of the set of allowed powers

        * Type 3: hash, 6 digit grid, power - ntype is negative.
        */

       if( (ntype >= 0) && (ntype <= 62) ) {
          int nu=ntype%10;
          if( nu == 0 || nu == 3 || nu == 7 ) {
             ndbm=ntype;
             memset(call_loc_pow,0,sizeof(char)*23);
             sprintf(cdbm,"%2d",ndbm);
             strncat(call_loc_pow,callsign,strlen(callsign));
             strncat(call_loc_pow," ",1);
             strncat(call_loc_pow,grid,4);
             strncat(call_loc_pow," ",1);
             strncat(call_loc_pow,cdbm,2);
             strncat(call_loc_pow,"\0",1);
             ihash=nhash(callsign,strlen(callsign),(uint32_t)146);
             strcpy(hashtab+ihash*13,callsign);
          } else {
             nadd=nu;
             if( nu > 3 ) nadd=nu-3;
             if( nu > 7 ) nadd=nu-7;
             n3=n2/128+32768*(nadd-1);
             if( !unpackpfx(n3,callsign) ) return 1;
             ndbm=ntype-nadd;
             memset(call_loc_pow,0,sizeof(char)*23);
             sprintf(cdbm,"%2d",ndbm);
             strncat(call_loc_pow,callsign,strlen(callsign));
             strncat(call_loc_pow," ",1);
             strncat(call_loc_pow,cdbm,2);
             strncat(call_loc_pow,"\0",1);
             int nu=ndbm%10;
             if( nu == 0 || nu == 3 || nu == 7 || nu == 10 ) {                //make sure power is OK
                ihash=nhash(callsign,strlen(callsign),(uint32_t)146);
                strcpy(hashtab+ihash*13,callsign);
             } else noprint=1;
          }
       } else if ( ntype < 0 ) {
          ndbm=-(ntype+1);
          memset(grid6,0,sizeof(char)*7);
    //        size_t len=strlen(callsign);
          size_t len=6;
          strncat(grid6,callsign+len-1,1);
          strncat(grid6,callsign,len-1);
          int nu=ndbm%10;
          if ((nu != 0 && nu != 3 && nu != 7 && nu != 10) ||
              !isalpha(grid6[0]) || !isalpha(grid6[1]) ||
              !isdigit(grid6[2]) || !isdigit(grid6[3])) {
             // not testing 4'th and 5'th chars because of this case: <PA0SKT/2> JO33 40
             // grid is only 4 chars even though this is a hashed callsign...
             //         isalpha(grid6[4]) && isalpha(grid6[5]) ) ) {
             noprint=1;
          }

          ihash=(n2-ntype-64)/128;
          if( strncmp(hashtab+ihash*13,"\0",1) != 0 ) {
             sprintf(callsign,"<%s>",hashtab+ihash*13);
          } else {
             sprintf(callsign,"%5s","<...>");
          }

          memset(call_loc_pow,0,sizeof(char)*23);
          sprintf(cdbm,"%2d",ndbm);
          strncat(call_loc_pow,callsign,strlen(callsign));
          strncat(call_loc_pow," ",1);
          strncat(call_loc_pow,grid6,strlen(grid6));
          strncat(call_loc_pow," ",1);
          strncat(call_loc_pow,cdbm,2);
          strncat(call_loc_pow,"\0",1);


          // I don't know what to do with these... They show up as "A000AA" grids.
          if( ntype == -64 ) noprint=1;
       }
       return noprint;
    }


   }   // namespace uwspr
} // namespace gr
