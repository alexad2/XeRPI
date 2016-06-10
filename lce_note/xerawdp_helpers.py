import ROOT
from collections import defaultdict
import pandas as pd


def load_xerawdp_tree(datasets,xerawdpPath):
    
    ROOT.gROOT.LoadMacro("stl_loader.h+")
    
    xerawdpTree = ROOT.TChain('T1')
    t2 = ROOT.TChain('T2')
    t3 = ROOT.TChain('T3')
    
    for dataset in datasets:
        xerawdpTree.AddFile(xerawdpPath+dataset+'.root')
        t2.AddFile(xerawdpPath+dataset+'.root')
        t3.AddFile(xerawdpPath+dataset+'.root')
    
    xerawdpTree.AddFriend(t2)
    xerawdpTree.AddFriend(t3)
    
    return xerawdpTree
    
def build_xerawdp_df(xerawdpTree):

    data = defaultdict(list)

    for i in range(xerawdpTree.GetEntries()):
        xerawdpTree.GetEntry(i)
        
        if xerawdpTree.NbS1Peaks < 2 or xerawdpTree.S1sTot[0] <= 0 or xerawdpTree.NbS2Peaks < 1:
            continue
    
        data['cs10Area'].append( xerawdpTree.cS1sTot[0] )
        data['s10Area'].append( xerawdpTree.S1sTot[0] )
        data['s10Time'].append( xerawdpTree.S1sPeak[0]*10.0 )
        data['s10Coin'].append( xerawdpTree.S1sCoin[0] )
        data['s10LeftEdge'].append( xerawdpTree.S1sLeftEdge[0]*10.0 )
        data['s10RightEdge'].append( xerawdpTree.S1sRightEdge[0]*10.0 )
    
        data['cs11Area'].append( xerawdpTree.cS1sTot[1] )
        data['s11Area'].append( xerawdpTree.S1sTot[1] )
        data['s11Time'].append( xerawdpTree.S1sPeak[1]*10.0 )
        data['s11Coin'].append( xerawdpTree.S1sCoin[1] )
        data['s11LeftEdge'].append( xerawdpTree.S1sLeftEdge[1]*10.0 )
        data['s11RightEdge'].append( xerawdpTree.S1sRightEdge[1]*10.0 )
        
        data['cs20Area'].append( xerawdpTree.cS2sTot[0] )
        data['s20Area'].append( xerawdpTree.S2sTot[0] )
        data['s20Time'].append( xerawdpTree.S2sPeak[0]*10.0 )
        data['s20Coin'].append( xerawdpTree.S2sCoin[0] )
        data['s20LeftEdge'].append( xerawdpTree.S2sLeftEdge[0]*10.0 )
        data['s20RightEdge'].append( xerawdpTree.S2sRightEdge[0]*10.0 )
        
        data['i0x'].append( xerawdpTree.cS2sPosNn[0][0]/10.0 )
        data['i0y'].append( xerawdpTree.cS2sPosNn[0][1]/10.0 )
        data['i0z'].append( -xerawdpTree.cS2sPosNn[0][2]/10.0 )
        
        try:
            data['cs21Area'].append( xerawdpTree.cS2sTot[1] )
            data['s21Area'].append( xerawdpTree.S2sTot[1] )
            data['s21Time'].append( xerawdpTree.S2sPeak[1]*10.0 )
            data['s21Coin'].append( xerawdpTree.S2sCoin[1] )
            data['s21LeftEdge'].append( xerawdpTree.S2sLeftEdge[1]*10.0 )
            data['s21RightEdge'].append( xerawdpTree.S2sRightEdge[1]*10.0 )
            
            data['i1x'].append( xerawdpTree.cS2sPosNn[1][0]/10.0 )
            data['i1y'].append( xerawdpTree.cS2sPosNn[1][1]/10.0 )
            data['i1z'].append( -xerawdpTree.cS2sPosNn[1][2]/10.0 )
        except IndexError:
            data['cs21Area'].append(0.0)
            data['s21Area'].append(0.0)
            data['s21Time'].append(0.0)
            data['s21Coin'].append(0.0)
            data['s21LeftEdge'].append( 0.0 )
            data['s21RightEdge'].append( 0.0 )
            
            data['i1x'].append( 15.0 )
            data['i1y'].append( 15.0 )
            data['i1z'].append( 35.0 )
            
    
    df_xerawdp = pd.DataFrame(data)
    return df_xerawdp
