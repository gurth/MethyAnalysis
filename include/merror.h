//
// Created by gurth on 4/6/21.
//

#ifndef METHYPROFILE_MERROR_H
#define METHYPROFILE_MERROR_H

#define SET_VALUE_ERROR 0x0F000
/* The highest bit of the value is 1 means numerical error.
 * Set value error.
 * Always happens when the format of bed file is not correct.
 * */

#define MUTEX_ERROR 0x00900
/* The highest bit of the second byte is 1 means system error.
 * Mutex init error.
 * */

#define ZLOG_ERROR 0x00a00
/* The highest bit of the second byte is 1 means system error.
 * ZLOG error.
 * */


#define FILE_OPEN_ERROR 0x00090
/* The highest bit of the third byte is 1 means I/O error.
 * File open, read, map, etc. error.
 * */

#define FILE_SAVE_ERROR 0x000A0
/* The highest bit of the third byte is 1 means I/O error.
 * File save error.
 * */

#define DIRECTORY_CREATE_ERROR 0x000d0
/* The highest bit of the third byte is 1 means I/O error.
 * Directory create error.
 * */

#define ARGUMENT_ERROR 0x00010
/* The highest bit of the third byte is 0 means common error.
 * Argument error.
 * */

#endif //METHYPROFILE_MERROR_H
