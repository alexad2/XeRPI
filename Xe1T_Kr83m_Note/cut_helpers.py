import ROOT
import numpy as  np

def cut_compare(df_xerawdp,df_pax,cut_info):
    
    key, cut_min, cut_max, bins, x_min, x_max, units = cut_info
    figure = './KrLce_Figures/f_'+key+'Hists.png'
    
    xerawdp_total = len(df_xerawdp.values)
    pax_total = len(df_pax.values)
    
    df_xerawdp_new = df_xerawdp
    df_pax_new = df_pax
    
    # Create and fill hists
    hists = []
    for i in range(2):
        hists.append(ROOT.TH1D('','',bins,x_min,x_max))
        if i==0:
            df = df_xerawdp
            title = 'Xerawdp'
        else:
            df = df_pax
            title = 'Pax'
        for j in range(len(df.values)):
            hists[i].Fill(df[key].values[j])
            
        hists[i].GetXaxis().SetTitle('%s (%s)'%(key,units))
        hists[i].GetXaxis().CenterTitle()
        hists[i].SetTitle(title+' '+key+' Histogram')
        hists[i].SetMinimum(10)
        hists[i].SetLineWidth(3)
            
    y_max = 1.2*max([hists[0].GetMaximum(),hists[1].GetMaximum()])
    
    cut_lines = []
    cut_str = ''
    if cut_min != 'none':
        
        cut_lines.append(ROOT.TLine(cut_min,0,cut_min,y_max))
        cut_lines[len(cut_lines)-1].SetLineColor(2)
        cut_lines[len(cut_lines)-1].SetLineWidth(3)
        
        cut_str += str(cut_min)+' <= '
        
        df_xerawdp_new = df_xerawdp_new[ df_xerawdp_new[key]>=cut_min ]
        df_pax_new = df_pax_new[ df_pax_new[key]>=cut_min ]
        
    cut_str += key
            
    if cut_max != 'none':
        
        cut_lines.append(ROOT.TLine(cut_max,0,cut_max,y_max))
        cut_lines[len(cut_lines)-1].SetLineColor(6)
        cut_lines[len(cut_lines)-1].SetLineWidth(3)
        
        cut_str += ' <= '+str(cut_max)
        
        df_xerawdp_new = df_xerawdp_new[ df_xerawdp_new[key]<=cut_max ]
        df_pax_new = df_pax_new[ df_pax_new[key]<=cut_max ]
                    
    xerawdp_acceptance = len(df_xerawdp_new.values)/xerawdp_total
    pax_acceptance = len(df_pax_new.values)/pax_total
    acceptance_ratios = [xerawdp_acceptance,pax_acceptance]
            
    c1 = ROOT.TCanvas('','',1600,700)
    ROOT.gStyle.SetOptStat(0)
    c1.Divide(2,1,0.02,0.02)
    
    pts = []
    
    for i in range(2):
        pad = c1.cd(i+1)
        pad.SetLogy()
        hists[i].SetMaximum(y_max)
        hists[i].Draw('hist')
        for line in cut_lines:
            line.Draw()
            
        pts.append(ROOT.TPaveText(.58, .68, .88, .88, 'NDC'))
        pts[i].AddText('Selection: '+cut_str)
        pts[i].AddText('Acceptance Ratio: %1.3f'%acceptance_ratios[i])
        pts[i].Draw()
        
    c1.Print(figure)
    c1.Clear()
    
    return figure

def apply_cuts(df,cuts):
    
    df_new = df
    
    for cut in cuts:
                    
        key, cut_min, cut_max, bins, x_min, x_max, units = cut
        if cut_min != 'none':
            df_new = df_new[ df_new[key]>=cut_min]
        if cut_max != 'none':
            df_new = df_new[ df_new[key]<=cut_max]
            
    return df_new

def n1_cuts(df,cuts):
    
    acceptances = np.zeros(len(cuts)+1)
    
    for i in range(len(cuts)+1):
        
        df_new = df
        
        for j, cut in enumerate(cuts):
            
            if i==j: continue
            
            key, cut_min, cut_max, bins, x_min, x_max, units = cut
            
            if cut_min != 'none':
                df_new = df_new[ df_new[key]>=cut_min ]
            if cut_max != 'none':
                df_new = df_new[ df_new[key]<=cut_max ]
                
        acceptances[i] = len(df_new.values)
        
    return acceptances
