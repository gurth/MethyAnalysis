//
// Created by gurth on 3/29/21.
//
// This file defines several macros which control the basic properties of this program.
// You can simply change their value to adapt particular situation, but make sure you
// follow the rules.
//

#ifndef METHYPROFILE_CONFIG_H
#define METHYPROFILE_CONFIG_H

//#define ALLOW_LOWERCASE

#define MAXTHREAD 16
/* Maximum thread number. Based on you computer performance. */

#define SEARCH_RANGE 512
/* Maximum search range, which avoid endless loop. This value should be larger than
 * maximum number of characters per line. */

#define SEARCH_RANGE_INDEX 102400  // 100 K
/* Maximum search range when search gene index.
 * Should be larger than maximum number of characters per entry in gff3 file. */

#define MAX_GENE 65535
/* Maximum gene number.
 * You can simply change the value if you sample more than 65536. */

#define MAX_CHR 128
/* Maximum chromosome number.
 * 128 fits most situations, but also we know a lot of species, such as
 * Acipenser sinensis (132), have more than 124 chromosomes. You should change this
 * value when this happens.
 * Notice that sex chromosome, mt chromosome and so on is coding as:
 *      X(Z)      Y(W)        MT         CH
 *   MAX_CHR-2  MAX_CHR-1  MAX_CHR-3  MAX_CHR-4
 * If there is also chloroplast chromosome data, please define macro __CHR_CHL.
 * If the sex chromosomes are Z and W, please define macro __CHR_ZW.
 * */

//#define __CHR_CHL    /* Define when chloroplast chromosome data is in sample. */
//#define __CHR_ZW     /* Define when sex chromosomes are Z and W. */

#define PROC_GENE_SIZE 0x80
/* Define how many genes can be processed per thread.
 * Can be larger if there are a large quantity of genes.*/

#define PROC_GENE_READ 0x800
/* Define how many genes can be processed at a time.
 * Please ensure: PROC_GENE_READ = PROC_GENE_SIZE x MAXTHREAD. */


//#define _FLAG_TEST
/* Define when you want to test the performance of this program.
 * Use this with Method::raw. */

#ifdef _FLAG_TEST

#define TEST_CHAR 'G'
/* Define in test. This program will fill the file with TEXT_CHAR.  */

#define BLOCK_SIZE 0x4000000 // 64 M
#define BLOCK_READ  0x40000000 // 1 G
/*Please see annotations below. */

#else

#define BLOCK_SIZE 0x400000 // 4 M
/* Define how many byte can be processed per thread. */

#define BLOCK_READ  0x4000000 // 64 M
/* Define how many byte can be processed at a time.
 * Please ensure: BLOCK_READ = BLOCK_SIZE x MAXTHREAD. */

#endif //!_FLAG_TEST

#define BLOCK_SIZE_INDEX 0x80000 // 512 K
/* Define how many byte can be processed per thread when search gene index. */

#define BLOCK_READ_INDEX 0x800000 // 4 M
/* Define how many byte can be processed at a time when search gene index.
 * Please ensure: BLOCK_READ_INDEX = BLOCK_SIZE_INDEX x MAXTHREAD. */

//#define SHOW_ALL_INFO
/* Show all information, which include how many bytes have been processed. */

#define SHOW_PROGRESSBAR
/* Show the progress bar. */

#define PROMOTER_LENGTH 2000
/* Promoter length at the front of gene. If you use -p option in arguments,
 * promoter methy fation will be output in profile.
 * Custom.
 * You can change it if needed. Default to 2000.
 * */

#define ALLOW_PLUG_IN_SAVE   /* Beta */
/* Define this if you want to use advanced feature of this program to save
 * profile. If you not sure about this, please ture it down.
 * */

#define CG_NUMBER
/* Define this when you want to save CG number in profile.
 * */

#define GENE_ID_BUFFER_LENGTH 0x20
/* Define gene ID buffer length.
 * Based on gene ID format.
 * Not important.
 * */

#define ENTRY_TYPE_BUFFER_LENGTH 0X20
/* Define the length of a entry type name. */

#define GO_BACK_LENGTH 0x10
/* When searching a string, the pointer may reach the edge of a block.
 * So the pointer need to go back.
 * The length must smaller than a line length.
 * */

#include "merror.h"

#endif //!METHYPROFILE_CONFIG_H
