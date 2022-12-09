import os
import simnibs.msh.transformations as transformations
import sys
import argparse
from simnibs import SIMNIBSDIR


def parseArguments(argv):
    # argument parsing exception handling
    parser = argparse.ArgumentParser(prog="electrode_warp",
                                     description= "maps eeg positions from MNI to subject space ")

    parser.add_argument("-c", metavar="EEG10-10_UI_Jurak_2007.csv", required=True,
                        help='Name of the eeg cap')
    parser.add_argument("-o", help="output base file name", metavar="10_10_positions", required=True)
    parser.add_argument("-s", help="subject directory (m2m folder)", metavar="m2m_ernie", required=True)

    args = parser.parse_args()
    if args.c is None:
        raise ValueError('a .csv-file has to be passed')

    if args.s is None:
        raise ValueError('subject directory needs to be passed')

    return args




def main():
    args = parseArguments(sys.argv)
    print("Cap Name:", args.c)
    print("Output Name:", args.o)
    print("Subject directory:", args.s)
    eeg_positions = os.path.join(args.s,"eeg_positions")
    if not os.path.exists(eeg_positions):
            os.mkdir(eeg_positions)
    cap_file = os.path.abspath(os.path.realpath(os.path.join(
                               SIMNIBSDIR,
                               'resources', 'ElectrodeCaps_MNI',
                               args.c)))
    cap_out = os.path.join(eeg_positions, args.o+'.csv')
    geo_out = os.path.join(eeg_positions, args.o+'.geo')
    transformations.warp_coordinates(
                    cap_file, args.s,
                    transformation_direction='mni2subject',
                    out_name=cap_out,
                    out_geo=geo_out)

if __name__ == '__main__':
    main()
