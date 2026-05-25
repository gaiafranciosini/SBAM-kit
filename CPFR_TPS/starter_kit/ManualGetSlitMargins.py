import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Example script that accepts width, height and an array of floats (exp)."
    )
    
    parser.add_argument(        
        "long",
        type=float,
        help="longer side of the slit [cm]"
    )

    parser.add_argument(        
        "short",
        type=float,
        help="shorter side of the slit [cm]"
    )


    parser.add_argument(
        "Wrx",
        type=float,
        help="Right (FLUKA) aperture width value [cm]"
    )

    parser.add_argument(
        "Wlx",
        type=float,
        help="Left (FLUKA) aperture width value [cm]"
    )

    parser.add_argument(
        "Htop",
        type=float,
        help="Top slit aperture height value [cm]"
    )

    parser.add_argument(
        "Hbot",
        type=float,
        help="Bottom aperture height value [cm]"
    )

    args = parser.parse_args()
    
    long=args.long
    short=args.short
    Wrx=args.Wrx
    Wlx=args.Wlx
    Htop=args.Htop
    Hbot=args.Hbot

#BEV x-axis and FLUKA x-axis are inverted

    Xi_right=Wrx
    Xf_right=Xi_right+short
    Xf_left=-Wlx
    Xi_left=Xf_left-short
    Yf_down=-Hbot
    Yi_down=Yf_down-short
    Yi_up=Htop
    Yf_up=Yi_up+short

    print("APERTURE: "+str(Xi_right)+" "+str(Xf_right)+" "+str(Xi_left)+" "+str(Xf_left)+" "+str(Yi_down)+" "+str(Yf_down)+" "+str(Yi_up)+" "+str(Yf_up))
    print("FIELD_SIZE: "+str(Wlx+Wrx)+" "+str(Htop+Hbot))
if __name__ == "__main__":
    main()
