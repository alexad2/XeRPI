import numpy as np
from collections import defaultdict
import ROOT
from subprocess import call
import pandas as pd

##################################################################################################

def atan(y, x):
    phi = np.arctan2(y, x)
    for i in range(len(phi)):
        if phi[i] < 0: 
            phi[i] += 2*np.pi
    return phi

##################################################################################################

def xe_to_lyBins(df,bin_settings,peak,bin_spec_dir='Bin_Hists'):
    
    R, Z, A_r, N_phi, N_z = bin_settings
    z_width = Z / N_z
    phi_widths = []
    for n in N_phi:
        phi_widths.append(2*np.pi/n)
    
    if peak == 's10':
        s1_spec_max = 20000
        s1_ene = 32.1498 # from nuclear data sheets a=83
        position = 'i0'
    elif peak == 's11':
        s1_spec_max = 100
        s1_ene = 9.4051
        position = 'i0' # only ever consider position of 1st s1
    else:
        print('error: invalid peak')
        return()

    bin_data = defaultdict(list)


    for z_i in range(int(N_z)):
    
        z_min = z_i * z_width
        z_max = (z_i+1) * z_width
        df_z = df[ (df[position+'z']<z_min) & (df[position+'z']>=z_max) ]
    
        for r_i in range(len(A_r)):
        
            if r_i == 0:
                r_min = 0
            else:
                r_min = A_r[r_i-1]
            r_max = A_r[r_i]
            df_r = df_z[ ( np.sqrt(df_z[position+'x']**2 + df_z[position+'y']**2)>r_min )
                       & ( np.sqrt(df_z[position+'x']**2 + df_z[position+'y']**2)<=r_max )]
        
            for phi_i in range(N_phi[r_i]):
            
                bin_data['z_i'].append(z_i)
                bin_data['z'].append( (z_max + z_min)/2 )
                bin_data['r_i'].append(r_i)
                bin_data['r'].append( (r_max + r_min)/2 )
                bin_data['phi_i'].append(phi_i)
                
                phi_min = phi_i * phi_widths[r_i] 
                phi_max = (phi_i+1) * phi_widths[r_i] 
                bin_data['phi'].append( (phi_max + phi_min)/2 )
            
                df_phi = df_r[ (atan(df_r[position+'y'].values, df_r[position+'x'].values) > phi_min) 
                             & (atan(df_r[position+'y'].values, df_r[position+'x'].values) <= phi_max )]
            
                bin_data['N'].append(len(df_phi))
                
                c1 = ROOT.TCanvas('','', 800, 700)
                hist = ROOT.TH1D('','', 100, 0, s1_spec_max)
                for i in range(len(df_phi[peak+'Area'])):
                    hist.Fill(df_phi[peak+'Area'].values[i])
                    
                if hist.GetEntries() < 1:
                    bin_data['ly'].append(-1)
                    bin_data['errly'].append(-1)
                    
                    bin_data['S1AreaMean'].append(-1)
                    bin_data['S1AreaMeanError'].append(-1)
                    
                    continue
                    
                hist.SetTitle(peak+' Spectrum: \
                              %.1f > z > %.1f, %.1f < r < %.1f, %.1f < phi < %.1f,'
                              %(z_min, z_max, r_min, r_max, phi_min, phi_max))
                hist.GetXaxis().SetTitle(peak+'Area (pe)')
                hist.GetXaxis().CenterTitle()
                hist.Sumw2()
                hist.SetStats(False)
                hist.Draw()
                hist.Fit('gaus')
                fit = hist.GetFunction('gaus')
                p1 = fit.GetParameter(1)
                e1 = fit.GetParError(1)
                
                bin_data['S1AreaMean'].append(p1)
                bin_data['S1AreaMeanError'].append(e1)
                
                bin_data['ly'].append(p1/s1_ene)
                bin_data['errly'].append(e1/s1_ene)
                
                if bin_spec_dir != 'none':
                    
                    call('mkdir '+bin_spec_dir,shell=True)
                
                    chi2 = fit.GetChisquare()
                    ndf = fit.GetNDF()
                    p0 = fit.GetParameter(0)
                    e0 = fit.GetParError(0)
                    p2 = fit.GetParameter(2)
                    e2 = fit.GetParError(2)
            
                    pt = ROOT.TPaveText(.58, .68, .88, .88, 'NDC')
                    pt.AddText('Entries = %d'%len(df_phi))
                    pt.AddText('#mu = %1.3f #pm %1.3f'%(p1, e1))
                    pt.AddText('#sigma = %1.3f #pm %1.3f' %(p2, e2))
                    pt.AddText('Amplitude = %1.3f #pm %1.3f' %(p0, e0))
                    pt.AddText('#chi^{2}/NDF = %1.3f/%1.3f' %(chi2, ndf))
            
                    pt.Draw()
            
                    c1.Print(bin_spec_dir+'/f_'+peak+'_z%d_r%d_phi%d.png' %(z_i, r_i, phi_i))
                
                c1.Clear()
                hist.Delete()
                
    return bin_data

##################################################################################################

def lyBins_to_txt(bin_data,out_file):
    
    f = open(out_file, 'w')
    header = 'z    t    r    zmid    tmid    rmid    ly    errly\n'
    f.write(header)
    
    for i in range(len(bin_data['z'])):
        bin_values = (str(bin_data['z_i'][i])+'    '+str(bin_data['phi_i'][i])+'    '+str(bin_data['r_i'][i])+'    '
                    +str(bin_data['z'][i])+'    '+str(bin_data['phi'][i])+'    '+str(bin_data['r'][i])+'    '
                    +str(bin_data['ly'][i])+'    '+str(bin_data['errly'][i])+'\n')
        f.write(bin_values)
                        
    f.close()    
    return
 
##################################################################################################

def bins_to_plot(bin_dict, peak, bin_settings, outfile, diff = False):
    if False:
        x=1
    else:     
        df = pd.DataFrame(bin_dict)
        dummyH_list=[]

        c1 = ROOT.TCanvas( '','', 2400, 3200 ) 
        ROOT.gStyle.SetOptStat(0)
        c1.Divide(3,4,0.02,0.02)

        z_hists = []
        
        max_ly = max(df['ly'])
        min_ly = min(df['ly'])
        
        zjump = bin_settings[1]/bin_settings[4]

        for z_i in range(int(bin_settings[4])):
            dummyH_list.append(ROOT.TH2D("","",100,-1*bin_settings[0],bin_settings[0],100,-1*bin_settings[0],bin_settings[0]))

            df_new = df[ df['z_i'] == z_i ]
            r_hists = []
            for r_i in range(len(bin_settings[2])):
                r_hists.append(ROOT.TH2D('','', bin_settings[3][r_i], 0, 2*np.pi, len(bin_settings[2]), 0, bin_settings[0]))
                df_newer = df_new[ df_new['r_i'] == r_i ]
                for i in range(len(df_newer)):
                    r_hists[r_i].Fill(df_newer['phi'].values[i], df_newer['r'].values[i],
                                    df_newer['ly'].values[i] )

            z_hists.append(r_hists)

            c1.cd(z_i+1)

            dummyH_list[z_i].Draw('colz')
            dummyH_list[z_i].SetTitle("%.2fcm < z < %.2fcm" %(z_i*zjump, (z_i+1)*zjump ))
            dummyH_list[z_i].GetZaxis().SetTitle("<s1Area>")
            dummyH_list[z_i].GetXaxis().SetTitle("x position [cm]")
            dummyH_list[z_i].GetXaxis().CenterTitle()
            dummyH_list[z_i].GetYaxis().SetTitle("y position [cm]")
            dummyH_list[z_i].GetYaxis().CenterTitle()

        #     c1.SetTopMargin(0.2)
            c1.SetRightMargin(0.2)

            for i in range(len(z_hists[z_i])):

                z_hists[z_i][i].GetZaxis().SetRangeUser(0, max_ly)
                if diff:
                    z_hists[z_i][i].GetZaxis().SetTitle("(pax_ly - xerawdp_ly)^{2} [pe/keV]")
                else:
                    z_hists[z_i][i].GetZaxis().SetTitle(peak + " ly [pe/keV]")
                z_hists[z_i][i].GetZaxis().SetTitleOffset(1.8)
                z_hists[z_i][i].Draw('pol colz a same') 

        c1.Print(outfile)
        c1.Clear()
        return

    ##################################################################################################

