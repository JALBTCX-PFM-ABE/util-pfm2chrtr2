
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.
*********************************************************************************************/


/****************************************  IMPORTANT NOTE  **********************************

    Comments in this file that start with / * ! are being used by Doxygen to document the
    software.  Dashes in these comment blocks are used to create bullet lists.  The lack of
    blank lines after a block of dash preceeded comments means that the next block of dash
    preceeded comments is a new, indented bullet list.  I've tried to keep the Doxygen
    formatting to a minimum but there are some other items (like <br> and <pre>) that need
    to be left alone.  If you see a comment that starts with / * ! and there is something
    that looks a bit weird it is probably due to some arcane Doxygen syntax.  Be very
    careful modifying blocks of Doxygen comments.

*****************************************  IMPORTANT NOTE  **********************************/




#ifndef VERSION

#define     VERSION     "PFM Software - pfm2chrtr2 V3.07 - 07/23/14"

#endif

/*

    Version 1.00
    Jan C. Depner
    08/12/10

    First working version.


    Version 1.01
    Jan C. Depner
    09/23/10

    Changed lower limit of total uncertainty to 0.0.  Apparently min_standard_dev isn't being populated properly
    in the PFM library...  DOH!!!


    Version 2.00
    Jan C. Depner
    11/22/10

    Added ability to run MISP or GMT to fill all null areas of the CHRTR2 file.


    Version 2.01
    Jan C. Depner
    02/23/11

    Cleaned up code just a bit.  Removed a couple of useless counters.


    Version 2.02
    Jan C. Depner
    05/06/11

    Fixed problem with getopt that only happens on Windows.


    Version 2.03
    Jan C. Depner
    11/18/11

    Removed GMT gridding option since no one was using (or wanted) it.


    Version 2.04
    Jan C. Depner
    01/24/12

    Correctly handle hand drawn contours.


    Version 2.05
    Jan C. Depner
    02/01/12

    Clear the CHRTR2_HEADER structure prior to populating it.


    Version 3.0
    Stacy Johnson
    06/20/2012

    Since chrtr2 was moved to grid registration we no longer need half node shift


    Version 3.01
    Stacy Johnson
    08/10/2012

    Modified to make chrtr2 creation based on the center of the pfm bin


    Version 3.02
    Stacy Johnson
    10/03/2012

    Made misp product run based on row/col as opposed to lat/lon as it is easier to track the process.
    Set FILTER to 0 as we want this conversion process to run over the entire file including edges.


    Version 3.03
    Stacy Johnson
    01/17/2013

    Multiply standard devation by 2 to keep inline with DUES 
    Changed code to always try to get uncertainty from the pfm if requested
    Changed code to change uncertainty of 0.0 to CHRTR2_NULL_Z_VALUE     


    Version 3.04
    Stacy Johnson
    02/08/2013

    Capped uncertainty as a user define percent of depth, defaults to 50%


    Version 3.05
    Stacy Johnson
    02/14/2013

    If H/V exceeds bounds set H/V to NULL in the chrtr2 for now until a better answer is determined
    Set depth to NULL if 0.0 


    Version 3.06
    Jan C. Depner
    03/08/14

    Added a line to force the linker to use -lm on Ubuntu 12.04 LTS x86_64 using
    gcc 4.6.3.  See bottom of main.c for complete explanation.


    Version 3.07
    Jan C. Depner (PFM Software)
    07/23/14

    - Switched from using the old NV_INT64 and NV_U_INT32 type definitions to the C99 standard stdint.h and
      inttypes.h sized data types (e.g. int64_t and uint32_t).

*/
