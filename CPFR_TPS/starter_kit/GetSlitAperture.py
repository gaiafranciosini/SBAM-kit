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
        "W",
        type=float,
        help="Aperture width value [cm]"
    )

    parser.add_argument(
        "H",
        type=float,
        help="Aperture height value [cm]"
    )

    parser.add_argument(
        "ape",
        type=float,
        help="Array of float values [mm] (ATTENTION mm !!!) "
    )

    args = parser.parse_args()
    
    long=args.long
    short=args.short
    W=args.W
    H=args.H
    ape=0.1*args.ape #mm->cm

    Xi_right=W*0.5+ape
    Xf_right=Xi_right+short
    Xf_left=-W*0.5-ape
    Xi_left=Xf_left-short
    Yf_down=-H*0.5-ape
    Yi_down=Yf_down-short
    Yi_up=H*0.5+ape
    Yf_up=Yi_up+short

    print("APERTURE: "+str(Xi_right)+" "+str(Xf_right)+" "+str(Xi_left)+" "+str(Xf_left)+" "+str(Yi_down)+" "+str(Yf_down)+" "+str(Yi_up)+" "+str(Yf_up))

if __name__ == "__main__":
    main()
