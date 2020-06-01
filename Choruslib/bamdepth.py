from Choruslib import Baminfo
from Choruslib import subprocesspath
import pysam


def continueregion(points, minlength=2):

    try:

        points.sort()

        start_index = 0

        end_index = 0

        continue_region = list()

        for index_now in range(1, len(points)):

            pre_index = index_now - 1

            if points[pre_index] + 1 == points[index_now]:

                if index_now == len(points) - 1:

                    if points[index_now] - points[start_index] + 1 >= minlength:
                        # print (points[start_index], points[index_now])
                        region_now = dict()
                        region_now['start_site'] = points[start_index]
                        region_now['end_site'] = points[index_now]
                        continue_region.append(region_now)

                else:

                    end_index = index_now

            else:

                if points[end_index] - points[start_index] + 1 >= minlength:
                    # print (points[start_index], points[end_index])
                    region_now = dict()
                    region_now['start_site'] = points[start_index]
                    region_now['end_site'] = points[end_index]
                    continue_region.append(region_now)

                start_index = index_now

                end_index = index_now

        return continue_region

    except Exception as e:

        print(('got exception in Jazzlib.region.continueregion: %r, terminating the pool' % (e,)))




def bamdepthtobed(bamfile, outbed='mindepth.bed', mindepth=1, minlength=1):

    baminfor = Baminfo.Baminfo(bamfile)

    outio = open(outbed,'w')

    for chrom in baminfor.getchrlen():

        print(chrom)

        depthstr = pysam.depth('-r', chrom, bamfile)

        depthchr = depthstr.split('\n')

        del depthstr

        points = list()

        for reg in depthchr:
            try:
                (chrom, site, depthnow) = reg.split('\t')
                site = int(site)
                depthnow = int(depthnow)

                if depthnow >= mindepth:
                    points.append(site)
            # print(chrom, site, depthnow)
            except Exception as e:
                print("warnning:", reg)



        continue_region = continueregion(points, minlength)

        for nowregion in continue_region:
            print(chrom, nowregion['start_site'],nowregion['end_site'], sep='\t', file=outio)

    outio.close()

# if __name__ == '__main__':
#
#     bamdepthtobed(bamfile='/mnt/e/Data/Solanum/Solanum_etuberosum/Solanum_etuberosum_map_to_DM.bam',
#                   outbed='/mnt/e/Data/Solanum/Solanum_etuberosum/Solanum_etuberosum_map_to_DM_min_depth.bed',
#                   mindepth=3,
#                   minlength=150)