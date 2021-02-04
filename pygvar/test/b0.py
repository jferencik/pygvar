from pygvar.reader import parse_gvar
import os
import pylab

if __name__ == '__main__':
    sample_gvars_folder = '/home/jano/Desktop/goes/sample_all_gvars/*'
    #sample_gvars_folder = '/media/d/data/goesw/OPERATIONAL/TO_PROCESS'
    src_folder = '/mnt/temp/disk2/2000/02/*goes08.2000.037.14*'
    from pygvar.navigation import le2ll

    import glob
    files = sorted(glob.glob(src_folder))
    print (files)

    for gvar_fpath in files:
        if not (gvar_fpath.endswith('.indx') or gvar_fpath.endswith('.meta')):
            _,gvar_name = os.path.split(gvar_fpath)

            chn = 4
            b0, ddict = parse_gvar(gvar_fpath, channels_to_extract=[chn])
            print (gvar_name, b0.scan_duration, b0.theoretical_scan_duration, b0.channels_shape.keys())
            sil, eil, sie, eie = b0.vis_bbox
            print (b0.igvln, b0.igvpx, b0.vis_bbox)#, [e for e in b0.coregistration_tbl_id], b0.index_of_active_corr_terms
            print (sil + (eil-sil)//2)


            #print [e for e in  b0.east_west_hourly_corr_terms]

            # g1cnt, g2cnt = b0.g1cnt, b0.g2cnt
            # grid1_detector = [e for e in b0.grid1_detector][:g1cnt]
            # grid2_detector = [e for e in b0.grid2_detector][:g2cnt]
            # grid1_pixel = [e for e in b0.grid1_pixel][:g1cnt]
            # grid2_pixel = [e for e in b0.grid2_pixel][:g2cnt]
            #
            #for e in zip(grid1_detector, grid1_pixel):


            d = ddict[chn]
            pylab.title(gvar_name)
            pylab.imshow(d)
            pylab.show()

