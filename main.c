
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



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <errno.h>
#include <string.h>
#include <getopt.h>

#include "nvutility.h"
#include "globals.hpp"

#include "pfm.h"
#include "chrtr2.h"
#include "chrtr2_shared.h"
#include "misp.h"

#include "version.h"

/*

    This program will take the average/MISP depth layer in a PFM file and create a CHRTR2 file.
    Horizontal and vertical uncertainty will be computed from an average of the horizontal and
    vertical uncertainty for all valid points in each bin.

    08/12/10
    Jan C. Depner

*/

#define         FILTER 0


void usage ()
{
  fprintf (stderr, "\nUsage: pfm2chrtr2 uncertainty_bound [--no_uncertainty] [--grid_type GRID_TYPE] [--output_file CHRTR2_FILE] PFM_FILE\n\n");
  fprintf (stderr, "\tWhere:\n\n");
  fprintf (stderr, "\t--no_uncertainty eliminates H/V uncertainty (but not total\n");
  fprintf (stderr, "\t\tuncertainty) from being stored in the output file.\n\n");
  fprintf (stderr, "\t--grid_type specifies the grid type where GRID_TYPE is\n");
  fprintf (stderr, "\t\tM or N for MISP, or NONE respectively.  If you\n");
  fprintf (stderr, "\t\tdo not specify the grid type the default is MISP.\n\n");
  fprintf (stderr, "\t--output_file specifies an output file name.  If you do\n");
  fprintf (stderr, "\t\tnot specify a name the output file will be the same as\n");
  fprintf (stderr, "\t\tthe PFM_FILE with the .pfm extension replaced with .ch2.\n");
  fprintf (stderr, "\tuncertainty_bound specifies the maximum uncertainty value\n");
  fprintf (stderr, "\t\tas a percentage of depth.\n\n\n");
  exit (-1);
}



static void add_point (NV_F64_COORD3 **xyz_array, NV_F64_COORD3 xyz, int32_t *count)
{
  *xyz_array = (NV_F64_COORD3 *) realloc (*xyz_array, (*count + 1) * sizeof (NV_F64_COORD3));
  if (*xyz_array == NULL)
    {
      perror ("Allocating xyz_array in remisp.");
      exit (-1);
    }


  (*xyz_array)[*count] = xyz;

  (*count)++;
}



/*  This function runs MISP on the selected area.  */

static void misp (int32_t weight, int32_t chrtr2_handle, CHRTR2_HEADER chrtr2_header)
{
  NV_F64_COORD3      *xyz_array = NULL, xyz;
  int32_t            i, j, out_count = 0, misp_weight;
  CHRTR2_RECORD      chrtr2_record;
  NV_F64_XYMBR       new_mbr;
  int32_t            gridcols, gridrows;
  float              *array = NULL;



  misp_weight = weight;


  /*  Number of rows and columns in the area  */

  gridcols = chrtr2_header.width;
  gridrows = chrtr2_header.height;


  /*  Save the data to memory.  */

  for (i = 0 ; i < gridrows ; i++)
    {
      for (j = 0 ; j < gridcols ; j++)
        {
	  chrtr2_read_record_row_col(chrtr2_handle,i,j,&chrtr2_record);


	  xyz.y = chrtr2_header.mbr.slat + (((float) i) * chrtr2_header.lat_grid_size_degrees);
	  xyz.x = chrtr2_header.mbr.wlon + (((float) j) * chrtr2_header.lon_grid_size_degrees);			
          xyz.z = chrtr2_record.z;


          /*  If we have data in the bin, go get it (we want to interpolate over already interpolated data  */
          /*  so we only load real or drawn data except in the filter border).  */

          if (chrtr2_record.status & (CHRTR2_REAL | CHRTR2_DIGITIZED_CONTOUR))
            {
              add_point (&xyz_array, xyz, &out_count);
            }
        }
    }


  /*  Don't process if we didn't have any input data (xyz would not have been allocated so we don't need to free it).  */

  if (!out_count)
    {
      fprintf (stderr, "\n\nNo data points found for gridding!\n\n");
      exit (-1);
    }


  /*  We're going to let MISP handle everything in zero based units of the bin size.  That is, we subtract off the  */
  /*  west lon from longitudes then divide by the grid size in the X direction.  We do the same with the latitude using  */
  /*  the south latitude.  This will give us values that range from 0.0 to gridcols in longitude and 0.0 to  */
  /*  gridrows in latitude.  */

  new_mbr.min_x = 0.0;
  new_mbr.min_y = 0.0;
  new_mbr.max_x = (double) gridcols;
  new_mbr.max_y = (double) gridrows;


  /*  Initialize the MISP engine.  */

  misp_init (1.0, 1.0, 0.05, 4, 20.0, 20, 999999.0, -999999.0, misp_weight, new_mbr);


  for (i = 0 ; i < out_count ; i++)
    {
      /*  Load the points.  */

      /*  IMPORTANT NOTE:  MISP (by default) grids using corner posts.  That is, the data in a bin is assigned to the 
          lower left corner of the bin.  Normal gridding/binning systems use the center of the bin.  Because of this we need
          to lie to MISP and tell it that the point is really half a bin lower and to the left.  This is extremely
          confusing but it works ;-)  */

      xyz.x = NINT((xyz_array[i].x - chrtr2_header.mbr.wlon) / chrtr2_header.lon_grid_size_degrees);
      xyz.y = NINT((xyz_array[i].y - chrtr2_header.mbr.slat) / chrtr2_header.lat_grid_size_degrees);
      xyz.z = xyz_array[i].z;
      misp_load (xyz);
    }


  fprintf (stderr, "Computing MISP surface\n");
  fflush (stderr);


  misp_proc ();

  fprintf (stderr, "Retrieving MISP data\n");
  fflush (stderr);


  /*  Allocating one more than gridcols due to constraints of old chrtr (see comments in misp_funcs.c)  */

  array = (float *) malloc ((gridcols + 1) * sizeof (float));

  if (array == NULL)
    {
      perror ("Allocating array in remisp");
      exit (-1);
    }


  /*  This is where we stuff the new interpolated surface back in to the CHRTR2.  */

  for (i = 0 ; i < gridrows ; i++)
    {
      if (!misp_rtrv (array)) break;


      for (j = 0 ; j < gridcols ; j++)
        {
          /*  Read the record.  */

	  chrtr2_read_record_row_col(chrtr2_handle,i,j,&chrtr2_record);


          /*  Only replace NULL values.  */

          if (chrtr2_record.status == CHRTR2_NULL)
            {
              /*  Mark the record as interpolated.  */

              chrtr2_record.status |= CHRTR2_INTERPOLATED;


              /*  If we exceeded the CHRTR2 limits we have to set it to the null depth (by definition, one greater than the max).  */

              if (array[j] <= chrtr2_header.max_z && array[j] >= chrtr2_header.min_z)
                {
                  chrtr2_record.z = array[j];
                }
              else
                {
                  chrtr2_record.z = chrtr2_header.max_z + 1.0;
                }


              /*  Write the record back out.  */

	      chrtr2_write_record_row_col(chrtr2_handle,i,j,chrtr2_record);
            }
        }
    }

  free (array);

  free (xyz_array);
}



int32_t main (int32_t argc, char *argv[])
{
  int32_t             i, j, k, numrecs, pfm_handle = 0, chrtr2_handle = 0, percent = 0, old_percent = -1, count, option_index, grid_type, ubound = 50;
  double              sum = 0.0, v_sum = 0.0, h_sum = 0.0;
  float               min_z, max_z;
  NV_I32_COORD2       coord;
  uint8_t             uncertainty;
  PFM_OPEN_ARGS       open_args;
  BIN_RECORD          bin_record;
  DEPTH_RECORD        *depth_record;
  CHRTR2_HEADER       chrtr2_header;
  CHRTR2_RECORD       chrtr2_record;
  char                c, chrtr2_file[512];
  extern char         *optarg;
  extern int          optind;


  fprintf(stderr, "\n\n %s \n\n", VERSION);
  fflush (stderr);


  option_index = 0;
  strcpy (chrtr2_file, "");
  uncertainty = NVTrue;
  grid_type = 1;

  while (NVTrue) 
    {
      static struct option long_options[] = {{"no_uncertainty", no_argument, 0, 0},
                                             {"grid_type", required_argument, 0, 0},
                                             {"output_file", required_argument, 0, 0},
                                             {0, no_argument, 0, 0}};

      c = (char) getopt_long (argc, argv, "", long_options, &option_index);
      if (c == -1) break;

      switch (c) 
        {
        case 0:

          switch (option_index)
            {
            case 0:
              uncertainty = NVFalse;
              break;

            case 1:
              if (strchr (optarg, 'N') || strchr (optarg, 'n'))
                {
                  grid_type = 0;
                }
              else if (strchr (optarg, 'M') || strchr (optarg, 'm'))
                {
                  grid_type = 1;
                }
              else if (strchr (optarg, 'G') || strchr (optarg, 'g'))
                {
                  grid_type = 2;
                }
              else
                {
                  usage ();
                }
              break;

            case 2:
              strcpy (chrtr2_file, optarg);
              break;

            case 3:
	      ubound = (100 / atoi (optarg));
	      break;
            }
          break;

        default:
          usage ();
          break;
        }
    }


  /*  Make sure we got the mandatory file name.  */

  if (optind >= argc) usage ();


  /*  Make sure it's the correct kind of file.  */

  if (!strstr (argv[optind], ".pfm")) usage ();


  strcpy (open_args.list_path, argv[optind]);


  /*  If the output file wasn't named on the command line, create it from the PFM file name.  */

  if (strlen (chrtr2_file) < 2)
    {
      strcpy (chrtr2_file, open_args.list_path);
      sprintf (&chrtr2_file[strlen (chrtr2_file) - 4], ".ch2");
    }
  else
    {
      /*  Make sure the .ch2 extension was included if the output file was specified on the command line.  */

      if (strcmp (&chrtr2_file[strlen (chrtr2_file) - 4], ".ch2")) strcat (chrtr2_file, ".ch2");
    }

  fprintf (stderr, "\n\nRejecting any uncertainty values greater than %d percent of depth\n\n", ubound);

  open_args.checkpoint = 0;
  pfm_handle = open_existing_pfm_file (&open_args);


  if (pfm_handle < 0) pfm_error_exit (pfm_error);


  /*  Check for projected data - yuck!  */

  if (open_args.head.proj_data.projection)
    {
      fprintf (stderr, "\n\npfm2chrtr2 does not handle projected data!!!\n\n");
      exit (-1);
    }


  /*  Populate the chrtr2 header prior to creating the file.  */

  memset (&chrtr2_header, 0, sizeof (CHRTR2_HEADER));

  strcpy (chrtr2_header.creation_software, VERSION);
  chrtr2_header.z_units = CHRTR2_METERS;
  chrtr2_header.mbr.wlon = open_args.head.mbr.min_x + (open_args.head.x_bin_size_degrees / 2);
  chrtr2_header.mbr.slat = open_args.head.mbr.min_y + (open_args.head.y_bin_size_degrees / 2);
  chrtr2_header.mbr.elon = open_args.head.mbr.max_x - (open_args.head.x_bin_size_degrees / 2);
  chrtr2_header.mbr.nlat = open_args.head.mbr.max_y - (open_args.head.y_bin_size_degrees / 2);

  chrtr2_header.width = open_args.head.bin_width;
  chrtr2_header.height = open_args.head.bin_height;
  chrtr2_header.lat_grid_size_degrees = open_args.head.y_bin_size_degrees;
  chrtr2_header.lon_grid_size_degrees = open_args.head.x_bin_size_degrees;
  chrtr2_header.min_z = -CHRTR2_NULL_Z_VALUE;
  chrtr2_header.max_z = CHRTR2_NULL_Z_VALUE;
  chrtr2_header.z_scale = open_args.scale;

  chrtr2_header.grid_type = CHRTR2_MISP;

  chrtr2_header.max_number_of_points = 16777215;
  chrtr2_header.min_uncertainty = 0.0;
  chrtr2_header.max_uncertainty = CHRTR2_NULL_Z_VALUE ; 
  chrtr2_header.uncertainty_scale = open_args.scale;
  strcpy (chrtr2_header.uncertainty_name, "Standard Deviation");

  if (uncertainty)
    {
      chrtr2_header.min_horizontal_uncertainty = 0.0;
      chrtr2_header.max_horizontal_uncertainty = 20000.0;
      chrtr2_header.horizontal_uncertainty_scale = open_args.head.horizontal_error_scale;
      chrtr2_header.min_vertical_uncertainty = 0.0;
      chrtr2_header.max_vertical_uncertainty = 10000.0;
      chrtr2_header.vertical_uncertainty_scale = open_args.head.vertical_error_scale;
    }
  else
    {
      chrtr2_header.horizontal_uncertainty_scale = 0.0;
      chrtr2_header.vertical_uncertainty_scale = 0.0;
    }


  /*  Try to create and open the chrtr2 file.  */

  chrtr2_handle = chrtr2_create_file (chrtr2_file, &chrtr2_header);
  if (chrtr2_handle < 0)
    {
      chrtr2_perror ();
      exit (-1);
    }


  min_z = 9999999999.0;
  max_z = -9999999999.0;


  /*  Loop through the PFM file.  */

  for (i = 0 ; i < open_args.head.bin_height ; i++)
    {
      coord.y = i;

      for (j = 0 ; j < open_args.head.bin_width ; j++)
        {
          coord.x = j;

          read_bin_record_index (pfm_handle, coord, &bin_record);

          memset (&chrtr2_record, 0, sizeof (CHRTR2_RECORD));

          if (bin_record.validity & PFM_DATA)
            {
              read_depth_array_index (pfm_handle, coord, &depth_record, &numrecs);

              sum = 0.0;
              v_sum = 0.0;
              h_sum = 0.0;
              count = 0;


              uint8_t drawn = NVFalse;
              for (k = 0 ; k < numrecs ; k++)
                {
                  if (!(depth_record[k].validity & (PFM_INVAL | PFM_DELETED | PFM_REFERENCE)))
                    {
                      //  Check for a hand-drawn contour (PFM_DATA is set in one or more of the depth records).

                      if (depth_record[k].validity & PFM_DATA) drawn = NVTrue;

                      if (uncertainty)
                        {
                          v_sum += depth_record[k].vertical_error;
                          h_sum += depth_record[k].horizontal_error;
                        }

                      sum += depth_record[k].xyz.z;
                      count++;
                    }
                }
              free (depth_record);


              /*  Just to be on the safe side let's make sure we got at least one valid point.  */

              if (count)
                {
                  if (uncertainty)
                    {
                      chrtr2_record.vertical_uncertainty = (float) (v_sum / (double) count);


                      /*  SJ - 02/12/2013 - temporarily set h to NULL when it exceeds the bounds  */

                      if (((float) (h_sum / (double) count)) >= chrtr2_header.max_horizontal_uncertainty)
			{
			  chrtr2_record.horizontal_uncertainty=chrtr2_header.max_horizontal_uncertainty-1;
			}
		      else
			{
			  chrtr2_record.horizontal_uncertainty = (float) (h_sum / (double) count);
			}
                    }

                  chrtr2_record.number_of_points = count;
                  chrtr2_record.uncertainty = bin_record.standard_dev * 2.0;
                  chrtr2_record.z = (sum / (double) count);


		  /*  SJ - 02/08/2013 - establish bound for uncertainty as a percent of depth  */

                  if ((bin_record.standard_dev * 2.0) > (chrtr2_record.z * ((float) ubound / 100.0))) 
		    {
                      chrtr2_record.uncertainty = CHRTR2_NULL_Z_VALUE;                         
		    }


		  //  SJ - 01/29/2013 0.0 is not a valid uncertainty.

                  if (chrtr2_record.uncertainty == 0.0) chrtr2_record.uncertainty = CHRTR2_NULL_Z_VALUE;


		  //  SJ - 02/14/2013 0.0 is not a valid depth, so set to NULL.

                  if (chrtr2_record.z == 0.0) chrtr2_record.z = CHRTR2_NULL_Z_VALUE;


                  if (drawn)
                    {
                      chrtr2_record.status = CHRTR2_DIGITIZED_CONTOUR;
                    }
                  else
                    {
                      chrtr2_record.status = CHRTR2_REAL;
                    }

                  min_z = MIN (chrtr2_record.z, min_z);
                  max_z = MAX (chrtr2_record.z, max_z);

                  if (chrtr2_write_record (chrtr2_handle, coord, chrtr2_record))
                    {
                      chrtr2_perror ();
                      exit (-1);
                    }
                }
            }
        }

      percent = ((float) i / (float) open_args.head.bin_height) * 100.0;
      if (percent != old_percent)
        {
          fprintf (stderr, "Processing - %03d%%\r", percent);
          fflush (stderr);
          old_percent = percent;
        }
    }

  printf("\n\n\n");

  chrtr2_header.min_observed_z = min_z;
  chrtr2_header.max_observed_z = max_z;

  chrtr2_update_header (chrtr2_handle, chrtr2_header);


  close_pfm_file (pfm_handle);


  chrtr2_close_file (chrtr2_handle);


  /*  MISP the data if requested.  */

  if (grid_type)
    {
      /*  Re-open the file and make sure it is a valid CHRTR2 file.  */

      chrtr2_handle = chrtr2_open_file (chrtr2_file, &chrtr2_header, CHRTR2_UPDATE);

      if (chrtr2_handle < 0)
        {
          fprintf (stderr, "The file %s is not a CHRTR2 structure or there was an error reading the file.\n", chrtr2_file);
          fprintf (stderr, "The error message returned was: %s\n\n", chrtr2_strerror ());

          exit (-1);
        }

      misp (2, chrtr2_handle, chrtr2_header);

      chrtr2_close_file (chrtr2_handle);
    }


  fprintf (stderr, "\nConversion complete\n\n");
  fflush (stderr);


  /*  Please ignore the following line.  It is useless.  Except...

      On some versions of Ubuntu, if I compile a program that doesn't use the math
      library but it calls a shared library that does use the math library I get undefined
      references to some of the math library functions even though I have -lm as the last
      library listed on the link line.  This happens whether I use qmake to build the
      Makefile or I have a pre-built Makefile.  Including math.h doesn't fix it either.
      The following line forces the linker to bring in the math library.  If there is a
      better solution please let me know at area.based.editor AT gmail DOT com.  */

  float ubuntu; ubuntu = 4.50 ; ubuntu = fmod (ubuntu, 1.0);


  return (0);
}
