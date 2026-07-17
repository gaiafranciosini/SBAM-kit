#!/usr/bin/env python3
import argparse
import os, sys
import numpy as np
import pydicom as dicom
from math import cos, sin, pi

# Configurazione Argomenti
parser = argparse.ArgumentParser(description='Export field info from DICOM RTPLAN')
parser.add_argument("file", help="dicom files", nargs='+', metavar='path')
args = parser.parse_args()

# Variabili globali
cts = []
rtstruct = []
rtplan = []
rtdose = []
beams = []

def loadDicoms():
    print('Caricamento file DICOM...')
    for fname in args.file:
        try:
            ds = dicom.dcmread(fname)
            uid = ds.SOPClassUID
            
            # CT
            if uid == '1.2.840.10008.5.1.4.1.1.2':
                cts.append(ds)
            # RT Structure Set
            elif uid == '1.2.840.10008.5.1.4.1.1.481.3':
                rtstruct.append(ds)
            # RT Plan (Fotoni) o RT Ion Plan (Protoni)
            elif uid in ['1.2.840.10008.5.1.4.1.1.481.5', '1.2.840.10008.5.1.4.1.1.481.8']:
                rtplan.append(ds)
            # RT Dose
            elif uid == '1.2.840.10008.5.1.4.1.1.481.2':
                rtdose.append(ds)
        except:
            continue

    print('CT trovate:', len(cts))
    print('RTPlan trovati:', len(rtplan))

def getBeamFrame(psideg, phideg, far):
    G0f = np.array([0,1,0])
    G0o = np.array([0,-far,0])
    G0u = np.array([0,0,1])
    
    psi = -pi/180 * psideg
    phi = pi/180 * phideg

    R_psi = np.array([
        [cos(psi),  sin(psi), 0],
        [-sin(psi), cos(psi), 0],
        [0,         0,        1]
    ])
    
    R_phi = np.array([
        [cos(phi),  0, sin(phi)],
        [0,         1, 0],
        [-sin(phi), 0, cos(phi)]
    ])

    Go = R_phi.dot(R_psi.dot(G0o))
    Gf = R_phi.dot(R_psi.dot(G0f))
    Gu = R_phi.dot(R_psi.dot(G0u))
    Gl = np.cross(Gu, Gf)

    return Go, Gf, Gu, Gl

def loadPLAN():
    if not rtplan:
        print('ERRORE: Nessun file RT Plan trovato.')
        return

    plan = rtplan[0]
    
    # Controlla se e un piano Ion (Protoni) o Standard
    is_ion = hasattr(plan, 'IonBeamSequence')
    beam_seq = plan.IonBeamSequence if is_ion else plan.BeamSequence
    
    print('Numero di fasci:', len(beam_seq))

    for beam in beam_seq:
        # Trova il primo Control Point
        cp = beam.IonControlPointSequence[0] if is_ion else beam.ControlPointSequence[0]
        
        gantry = float(getattr(cp, 'GantryAngle', 0))
        couch = float(getattr(cp, 'PatientSupportAngle', 0))
        iso = np.array(getattr(cp, 'IsocenterPosition', [0,0,0])) / 10.0
        
        # SAD (Source Axis Distance)
        sad = float(getattr(beam, 'SourceAxisDistance', 1000)) / 10.0
        
        Go, Gf, Gu, Gl = getBeamFrame(gantry, couch, sad)
        Go += iso

        # Pulizia zeri
        for v in [Go, Gf, Gu, Gl]:
            v[np.abs(v) < 1e-5] = 0.0

        beams.append({
            'ID': beam.BeamNumber,
            'O': Go, 'f': Gf, 'u': Gu, 'l': Gl, 'SAD': sad
        })

# Esecuzione
loadDicoms()
loadPLAN()

if beams:
    with open('beams.txt', 'w') as f:
        for b in beams:
            line = str(b['ID']) + " "
            vals = [b['O'][0], b['O'][1], b['O'][2], 
                    b['l'][0], b['l'][1], b['l'][2], 
                    b['u'][0], b['u'][1], b['u'][2], 
                    b['f'][0], b['f'][1], b['f'][2], 
                    b['SAD']]
            line += " ".join(["%.4f" % v for v in vals])
            print(line)
            f.write(line + "\n")
    print('File beams.txt generato con successo.')


