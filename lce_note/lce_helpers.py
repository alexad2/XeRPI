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

def xe100_to_lyBins(df,bin_settings,peak,bin_spec_dir='Bin_Hists'):
    
    R, Z, A_r, N_phi, N_z = bin_settings
    z_width = Z / N_z
    phi_widths = []
    for n in N_phi:
        phi_widths.append(2*np.pi/n)
    
    if peak == 's10':
        s1_spec_max = 200
        s1_ene = 32.1498 # from nuclear data sheets a=83
        postion = 'i0'
    elif peak == 's11':
        s1_spec_max = 100
        s1_ene = 9.4051
        position = 'i1'
    else:
        print('error: invalid peak')
        return()

    bin_data = defaultdict(list)


    for z_i in range(int(N_z)):
    
        z_min = z_i * z_width
        z_max = (z_i+1) * z_width
        df_z = df[ (df[position+'z']>z_min) & (df[position+'z']<=z_max) ]
    
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
                hist.SetTitle(peak+' Spectrum: \
                              %.1f < z < %.1f, %.1f < r < %.1f, %.1f < phi < %.1f,'
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

def txt_to_plot(file_name, peak, bin_settings, diff = False):
    if not isinstance(file_name, str): 
        print ('Error: arg must be str')
    else:
        fl = open(file_name, 'r')
        ls = fl.readline().strip().split(" ")
        labels = []
        for item in ls:
            if item != '':
                labels.append(item)
        fl.close()
        bin_array = np.loadtxt(file_name, skiprows=1)
        bin_array_new = bin_array.tolist()
        bin_array_new.sort(key = lambda x: (float(x[0]), float(x[2]), float(x[1])))
       
        for j in range(len(bin_array_new)):
            
            if bin_array_new[j][3] < -10:
                bin_array_new[j][3] = bin_array_new[j][3] / 10
                bin_array_new[j][5] = bin_array_new[j][5] / 10
            
            if bin_array_new[j][3] < 0:
                bin_array_new[j][3] = bin_array_new[j][3] * -1

        bin_dict = defaultdict(list)
        for i in range(len(bin_array_new[0])):
            for j in range(len(bin_array_new)):
                bin_dict[labels[i]].append(bin_array_new[j][i])
                
                
                
        df = pd.DataFrame(bin_dict)
        dummyH_list=[]

        c1 = ROOT.TCanvas( '','', 2400, 3200 ) 
        ROOT.gStyle.SetOptStat(0)
        c1.Divide(3,4,0.02,0.02)

        z_hists = []
        
        max_ly = max(df['ly'])
        min_ly = min(df['ly'])

        for z_i in range(int(bin_settings[4])):
            dummyH_list.append(ROOT.TH2D("","",100,-1*bin_settings[0],bin_settings[0],100,-1*bin_settings[0],bin_settings[0]))

            df_new = df[ df['z'] == z_i ]
            r_hists = []
            for r_i in range(len(bin_settings[2])):
                r_hists.append(ROOT.TH2D('','', bin_settings[3][r_i], 0, 2*np.pi, len(bin_settings[2]), 0, bin_settings[0]))
                df_newer = df_new[ df_new['r'] == r_i ]
                for i in range(len(df_newer)):
                    r_hists[r_i].Fill(df_newer['tmid'].values[i], df_newer['rmid'].values[i],
                                    df_newer['ly'].values[i] )

            z_hists.append(r_hists)

            c1.cd(z_i+1)

            dummyH_list[z_i].Draw('colz')
            dummyH_list[z_i].SetTitle("%.2fcm < z < %.2fcm" %(z_i*3.03, (z_i+1)*3.03 ))
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

        c1.Print(file_name.split('.')[0] + '_' + peak + '.png')
        c1.Clear()
        return
    
##################################################################################################

def difference_in_ly(file1, file2, outfile):
    
    # GIVES ABSOLUTE VALUE OF DIFFERENCES
    
    if not isinstance(file1, str) or not isinstance(file2, str): 
        print ('Error: args must be str type')
        
    else:
        
        # read first file
        
        fl = open(file1, 'r')
        ls = fl.readline().strip().split(" ")
        labels = []
        for item in ls:
            if item != '':
                labels.append(item)
        fl.close()
        bin_array = np.loadtxt(file1, skiprows=1)
        bin_array_new_1 = bin_array.tolist()
        bin_array_new_1.sort(key = lambda x: (float(x[0]), float(x[2]), float(x[1])))
       
        for j in range(len(bin_array_new_1)):
            
            if bin_array_new_1[j][3] < -10:
                bin_array_new_1[j][3] = bin_array_new_1[j][3] / 10
                bin_array_new_1[j][5] = bin_array_new_1[j][5] / 10
            
            if bin_array_new_1[j][3] < 0:
                bin_array_new_1[j][3] = bin_array_new_1[j][3] * -1

        bin_dict_1 = defaultdict(list)
        for i in range(len(bin_array_new_1[0])):
            for j in range(len(bin_array_new_1)):
                bin_dict_1[labels[i]].append(bin_array_new_1[j][i])
                
        # read second file
        
        fl = open(file2, 'r')
        ls = fl.readline().strip().split(" ")
        labels = []
        for item in ls:
            if item != '':
                labels.append(item)
        fl.close()
        bin_array = np.loadtxt(file2, skiprows=1)
        bin_array_new_2 = bin_array.tolist()
        bin_array_new_2.sort(key = lambda x: (float(x[0]), float(x[2]), float(x[1])))
       
        for j in range(len(bin_array_new_2)):
            
            if bin_array_new_2[j][3] < -10:
                bin_array_new_2[j][3] = bin_array_new_2[j][3] / 10
                bin_array_new_2[j][5] = bin_array_new_2[j][5] / 10
            
            if bin_array_new_2[j][3] < 0:
                bin_array_new_2[j][3] = bin_array_new_2[j][3] * -1

        bin_dict_2 = defaultdict(list)
        for i in range(len(bin_array_new_2[0])):
            for j in range(len(bin_array_new_2)):
                bin_dict_2[labels[i]].append(bin_array_new_2[j][i])
    
        
        
        f = open(outfile, 'w')
        f.write('z    t    r    zmid    tmid    rmid    ly    errly')
        for i in range(len(bin_dict_1['z'])):
            f.write('\n' + str(bin_dict_1['z'][i])+'    '+str(bin_dict_1['t'][i])+'    '+str(bin_dict_1['r'][i])+'    '
                        +str(bin_dict_1['zmid'][i])+'    '+str(bin_dict_1['tmid'][i])+'    '+str(bin_dict_1['rmid'][i])+'    '
                        +str(np.fabs(bin_dict_1['ly'][i]-bin_dict_2['ly'][i])) + '    '+'0.0')
        f.close()    
        return
    
##################################################################################################

# File 1: Xerawdp, File 2: Pax

def triple_plot(file1, file2, fileDiff, peak, bin_settings, figure_name):
    
    fl = open(file1, 'r')
    ls = fl.readline().strip().split(" ")
    labels = []
    for item in ls:
        if item != '':
            labels.append(item)
    fl.close()
    bin_array = np.loadtxt(file1, skiprows=1)
    bin_array_new = bin_array.tolist()
    bin_array_new.sort(key = lambda x: (float(x[0]), float(x[2]), float(x[1])))

    for j in range(len(bin_array_new)):

        if bin_array_new[j][3] < -10:
            bin_array_new[j][3] = bin_array_new[j][3] / 10
            bin_array_new[j][5] = bin_array_new[j][5] / 10

        if bin_array_new[j][3] < 0:
            bin_array_new[j][3] = bin_array_new[j][3] * -1

    bin_dict = defaultdict(list)
    for i in range(len(bin_array_new[0])):
        for j in range(len(bin_array_new)):
            bin_dict[labels[i]].append(bin_array_new[j][i])
    ###        
    df1 = pd.DataFrame(bin_dict)
    ###
    
    fl = open(file2, 'r')
    ls = fl.readline().strip().split(" ")
    labels = []
    for item in ls:
        if item != '':
            labels.append(item)
    fl.close()
    bin_array = np.loadtxt(file2, skiprows=1)
    bin_array_new = bin_array.tolist()
    bin_array_new.sort(key = lambda x: (float(x[0]), float(x[2]), float(x[1])))

    for j in range(len(bin_array_new)):

        if bin_array_new[j][3] < -10:
            bin_array_new[j][3] = bin_array_new[j][3] / 10
            bin_array_new[j][5] = bin_array_new[j][5] / 10

        if bin_array_new[j][3] < 0:
            bin_array_new[j][3] = bin_array_new[j][3] * -1

    bin_dict = defaultdict(list)
    for i in range(len(bin_array_new[0])):
        for j in range(len(bin_array_new)):
            bin_dict[labels[i]].append(bin_array_new[j][i])
    ###
    df2 = pd.DataFrame(bin_dict)
    ##
    
    fl = open(fileDiff, 'r')
    ls = fl.readline().strip().split(" ")
    labels = []
    for item in ls:
        if item != '':
            labels.append(item)
    fl.close()
    bin_array = np.loadtxt(fileDiff, skiprows=1)
    bin_array_new = bin_array.tolist()
    bin_array_new.sort(key = lambda x: (float(x[0]), float(x[2]), float(x[1])))

    for j in range(len(bin_array_new)):

        if bin_array_new[j][3] < -10:
            bin_array_new[j][3] = bin_array_new[j][3] / 10
            bin_array_new[j][5] = bin_array_new[j][5] / 10

        if bin_array_new[j][3] < 0:
            bin_array_new[j][3] = bin_array_new[j][3] * -1

    bin_dict = defaultdict(list)
    for i in range(len(bin_array_new[0])):
        for j in range(len(bin_array_new)):
            bin_dict[labels[i]].append(bin_array_new[j][i])
    ###
    dfDiff = pd.DataFrame(bin_dict)
    ##
    
    
    dummyH_list=[[],[],[]]
    
    c1 = ROOT.TCanvas( '','', 2400, 3200 ) 
    ROOT.gStyle.SetOptStat(0)
    c1.RangeAxis(0.0, 0.0, 12.0, 12.0)
    c1.Divide(6,8,0.02,0.02)
    ROOT.gStyle.SetTitleFontSize(0.1)
   
    hline1 = ROOT.TLine(0.0, 0.0, 0.33, 0.0)
    hline1.Draw()
    hline2 = ROOT.TLine(0.0, 0.25, 1.0, 0.25)
    hline2.Draw()
    hline3 = ROOT.TLine(0.0, 0.5, 1.0, 0.5)
    hline3.Draw()
    hline4 = ROOT.TLine(0.0, 0.75, 1.0, 0.75)
    hline4.Draw()
    hline5 = ROOT.TLine(0.0, 1.0, 1.0, 1.0)
    hline5.Draw()
    
    vline1 = ROOT.TLine(0.0, 0.0, 0.0, 1.0)
    vline1.Draw()
    vline2 = ROOT.TLine(0.33, 0.0, 0.33, 1.0)
    vline2.Draw()
    vline3 = ROOT.TLine(0.66, 0.25, 0.66, 1.0)
    vline3.Draw()
    vline4 = ROOT.TLine(1.0, 0.25, 1.0, 1.0)
    vline4.Draw()
    
    title1 = ROOT.TPaveText(0.11, 0.98, 0.22, 0.995)
    title1.AddText('0.00 cm < z < 3.03 cm')
    title1.Draw()
    title2 = ROOT.TPaveText(0.44, 0.98, 0.55, 0.995)
    title2.AddText('3.03 cm < z < 6.06 cm')
    title2.Draw()
    title3 = ROOT.TPaveText(0.77, 0.98, 0.88, 0.995)
    title3.AddText('6.06 cm < z < 9.09 cm')
    title3.Draw()
    title4 = ROOT.TPaveText(0.11, 0.73, 0.22, 0.745)
    title4.AddText('9.09 cm < z < 12.12 cm')
    title4.Draw()
    title5 = ROOT.TPaveText(0.44, 0.73, 0.55, 0.745)
    title5.AddText('12.12 cm < z < 15.15 cm')
    title5.Draw()
    title6 = ROOT.TPaveText(0.77, 0.73, 0.88, 0.745)
    title6.AddText('15.15 cm < z < 18.18 cm')
    title6.Draw()
    title7 = ROOT.TPaveText(0.11, 0.48, 0.22, 0.495)
    title7.AddText('18.18 cm < z < 21.21 cm')
    title7.Draw()
    title8 = ROOT.TPaveText(0.44, 0.48, 0.55, 0.495)
    title8.AddText('21.21 cm < z < 24.24 cm')
    title8.Draw()
    title9 = ROOT.TPaveText(0.77, 0.48, 0.88, 0.495)
    title9.AddText('24.24 cm < z < 27.27 cm')
    title9.Draw()
    title10 = ROOT.TPaveText(0.11, 0.23, 0.22, 0.245)
    title10.AddText('27.27 cm < z < 30.30 cm')
    title10.Draw()
    
    z_hists_1 = []
    z_hists_2 = []
    z_hists_diff = []
    
    for z_i in range(int(bin_settings[4])):
        
        dummyH_list[0].append(ROOT.TH2D("","",100,-1*bin_settings[0],bin_settings[0],100,-1*bin_settings[0],bin_settings[0]))
        dummyH_list[1].append(ROOT.TH2D("","",100,-1*bin_settings[0],bin_settings[0],100,-1*bin_settings[0],bin_settings[0]))
        dummyH_list[2].append(ROOT.TH2D("","",100,-1*bin_settings[0],bin_settings[0],100,-1*bin_settings[0],bin_settings[0]))

        df_new_1 = df1[ df1['z'] == z_i ]
        r_hists_1 = []
        df_new_2 = df2[ df2['z'] == z_i ]
        r_hists_2 = []
        df_new_diff = dfDiff[ dfDiff['z'] == z_i ]
        r_hists_diff = []
        for r_i in range(len(bin_settings[2])):

            r_hists_1.append(ROOT.TH2D('','', bin_settings[3][r_i], 0, 2*np.pi, len(bin_settings[2]), 0, bin_settings[0]))
            df_newer_1 = df_new_1[ df_new_1['r'] == r_i ]
            for i in range(len(df_newer_1)):
                r_hists_1[r_i].Fill(df_newer_1['tmid'].values[i], df_newer_1['rmid'].values[i],
                                df_newer_1['ly'].values[i] )

            r_hists_2.append(ROOT.TH2D('','', bin_settings[3][r_i], 0, 2*np.pi, len(bin_settings[2]), 0, bin_settings[0]))
            df_newer_2 = df_new_2[ df_new_2['r'] == r_i ]
            for i in range(len(df_newer_2)):
                r_hists_2[r_i].Fill(df_newer_2['tmid'].values[i], df_newer_2['rmid'].values[i],
                                df_newer_2['ly'].values[i] )

            r_hists_diff.append(ROOT.TH2D('','', bin_settings[3][r_i], 0, 2*np.pi, len(bin_settings[2]), 0, bin_settings[0]))
            df_newer_diff = df_new_diff[ df_new_diff['r'] == r_i ]
            for i in range(len(df_newer_diff)):
                r_hists_diff[r_i].Fill(df_newer_diff['tmid'].values[i], df_newer_diff['rmid'].values[i],
                                df_newer_diff['ly'].values[i] )


        z_hists_1.append(r_hists_1)
        z_hists_2.append(r_hists_2)
        z_hists_diff.append(r_hists_diff)

        if z_i < 3:
            padno = 2*z_i + 1
        elif z_i < 6:
            padno = 2*z_i + 7
        elif z_i < 9:
            padno = 2*z_i + 13
        elif z_i == 9:
            padno = 2*z_i + 19

        pad1 = c1.cd(padno)
        pad1.SetBorderSize(20)
        pad1.SetBorderMode(1)

        dummyH_list[0][z_i].Draw('colz')
        dummyH_list[0][z_i].SetTitle("Xerawdp") #%.2fcm < z < %.2fcm" %(z_i*3.03, (z_i+1)*3.03 
        dummyH_list[0][z_i].GetZaxis().SetTitle("<s1Area>")
        dummyH_list[0][z_i].GetXaxis().SetTitle("x position [cm]")
        dummyH_list[0][z_i].GetXaxis().SetTitleSize(0.06)
        dummyH_list[0][z_i].GetXaxis().CenterTitle()
        dummyH_list[0][z_i].GetYaxis().SetTitle("y position [cm]")
        dummyH_list[0][z_i].GetYaxis().SetTitleSize(0.06)
        dummyH_list[0][z_i].GetYaxis().CenterTitle()
        c1.SetRightMargin(0.2)

        for i in range(len(z_hists_1[z_i])):

            z_hists_1[z_i][i].GetZaxis().SetRangeUser(0, 7)
            z_hists_1[z_i][i].GetZaxis().SetTitle(peak + " ly [pe/keV]")
            z_hists_1[z_i][i].GetZaxis().SetTitleOffset(1.8)
            z_hists_1[z_i][i].Draw('pol colz a same')

        pad2 = c1.cd(padno+1)
        pad2.SetBorderSize(20)
        pad2.SetBorderMode(1)

        dummyH_list[1][z_i].Draw('colz')
        dummyH_list[1][z_i].SetTitle("Pax")
        dummyH_list[1][z_i].GetZaxis().SetTitle("<s1Area>")
        dummyH_list[1][z_i].GetXaxis().SetTitle("x position [cm]")
        dummyH_list[1][z_i].GetXaxis().SetTitleSize(0.06)
        dummyH_list[1][z_i].GetXaxis().CenterTitle()
        dummyH_list[1][z_i].GetYaxis().SetTitle("y position [cm]")
        dummyH_list[1][z_i].GetYaxis().SetTitleSize(0.06)
        dummyH_list[1][z_i].GetYaxis().CenterTitle()
        c1.SetRightMargin(0.2)

        for i in range(len(z_hists_2[z_i])):

            z_hists_2[z_i][i].GetZaxis().SetRangeUser(0, 7)
            z_hists_2[z_i][i].GetZaxis().SetTitle(peak + " ly [pe/keV]")
            z_hists_2[z_i][i].GetZaxis().SetTitleOffset(1.8)
            z_hists_2[z_i][i].Draw('pol colz a same')

        pad3 = c1.cd(padno+6)
        pad3.SetBorderSize(20)
        pad3.SetBorderMode(1)


        dummyH_list[2][z_i].Draw('colz')
        dummyH_list[2][z_i].SetTitle("|Difference|")
        dummyH_list[2][z_i].GetZaxis().SetTitle("<s1Area>")
        dummyH_list[2][z_i].GetXaxis().SetTitle("x position [cm]")
        dummyH_list[2][z_i].GetXaxis().SetTitleSize(0.06)
        dummyH_list[2][z_i].GetXaxis().CenterTitle()
        dummyH_list[2][z_i].GetYaxis().SetTitle("y position [cm]")
        dummyH_list[2][z_i].GetYaxis().SetTitleSize(0.06)
        dummyH_list[2][z_i].GetYaxis().CenterTitle()
        c1.SetRightMargin(0.2)

        for i in range(len(z_hists_diff[z_i])):

            z_hists_diff[z_i][i].GetZaxis().SetRangeUser(0, 1)
            z_hists_diff[z_i][i].GetZaxis().SetTitle(peak + " ly [pe/keV]")
            z_hists_diff[z_i][i].GetZaxis().SetTitleOffset(1.8)
            z_hists_diff[z_i][i].Draw('pol colz a same')
    
    for i in range(1,49):
        pad = c1.cd(i)
        pad.SetTopMargin(0.2)
        pad.SetBottomMargin(0)
    
    c1.Print(figure_name)
    c1.Clear()
    return

##################################################################################################

# df1: Xerawdp, df2: Pax

def df_to_plot (df1, df2, dfDiff, peak, bin_settings, figure_name):

    dummyH_list=[[],[],[]]
    
    c1 = ROOT.TCanvas( '','', 2400, 3200 ) 
    ROOT.gStyle.SetOptStat(0)
    c1.RangeAxis(0.0, 0.0, 12.0, 12.0)
    c1.Divide(6,8,0.02,0.02)
    ROOT.gStyle.SetTitleFontSize(0.1)
   
    hline1 = ROOT.TLine(0.0, 0.0, 0.33, 0.0)
    hline1.Draw()
    hline2 = ROOT.TLine(0.0, 0.25, 1.0, 0.25)
    hline2.Draw()
    hline3 = ROOT.TLine(0.0, 0.5, 1.0, 0.5)
    hline3.Draw()
    hline4 = ROOT.TLine(0.0, 0.75, 1.0, 0.75)
    hline4.Draw()
    hline5 = ROOT.TLine(0.0, 1.0, 1.0, 1.0)
    hline5.Draw()
    
    vline1 = ROOT.TLine(0.0, 0.0, 0.0, 1.0)
    vline1.Draw()
    vline2 = ROOT.TLine(0.33, 0.0, 0.33, 1.0)
    vline2.Draw()
    vline3 = ROOT.TLine(0.66, 0.25, 0.66, 1.0)
    vline3.Draw()
    vline4 = ROOT.TLine(1.0, 0.25, 1.0, 1.0)
    vline4.Draw()
    
    title1 = ROOT.TPaveText(0.11, 0.98, 0.22, 0.995)
    title1.AddText('0.00 cm < z < 3.03 cm')
    title1.Draw()
    title2 = ROOT.TPaveText(0.44, 0.98, 0.55, 0.995)
    title2.AddText('3.03 cm < z < 6.06 cm')
    title2.Draw()
    title3 = ROOT.TPaveText(0.77, 0.98, 0.88, 0.995)
    title3.AddText('6.06 cm < z < 9.09 cm')
    title3.Draw()
    title4 = ROOT.TPaveText(0.11, 0.73, 0.22, 0.745)
    title4.AddText('9.09 cm < z < 12.12 cm')
    title4.Draw()
    title5 = ROOT.TPaveText(0.44, 0.73, 0.55, 0.745)
    title5.AddText('12.12 cm < z < 15.15 cm')
    title5.Draw()
    title6 = ROOT.TPaveText(0.77, 0.73, 0.88, 0.745)
    title6.AddText('15.15 cm < z < 18.18 cm')
    title6.Draw()
    title7 = ROOT.TPaveText(0.11, 0.48, 0.22, 0.495)
    title7.AddText('18.18 cm < z < 21.21 cm')
    title7.Draw()
    title8 = ROOT.TPaveText(0.44, 0.48, 0.55, 0.495)
    title8.AddText('21.21 cm < z < 24.24 cm')
    title8.Draw()
    title9 = ROOT.TPaveText(0.77, 0.48, 0.88, 0.495)
    title9.AddText('24.24 cm < z < 27.27 cm')
    title9.Draw()
    title10 = ROOT.TPaveText(0.11, 0.23, 0.22, 0.245)
    title10.AddText('27.27 cm < z < 30.30 cm')
    title10.Draw()
    
    z_hists_1 = []
    z_hists_2 = []
    z_hists_diff = []
    
    for z_i in range(int(bin_settings[4])):
        
        dummyH_list[0].append(ROOT.TH2D("","",100,-1*bin_settings[0],bin_settings[0],100,-1*bin_settings[0],bin_settings[0]))
        dummyH_list[1].append(ROOT.TH2D("","",100,-1*bin_settings[0],bin_settings[0],100,-1*bin_settings[0],bin_settings[0]))
        dummyH_list[2].append(ROOT.TH2D("","",100,-1*bin_settings[0],bin_settings[0],100,-1*bin_settings[0],bin_settings[0]))

        df_new_1 = df1[ df1['z'] == z_i ]
        r_hists_1 = []
        df_new_2 = df2[ df2['z'] == z_i ]
        r_hists_2 = []
        df_new_diff = dfDiff[ dfDiff['z'] == z_i ]
        r_hists_diff = []
        for r_i in range(len(bin_settings[2])):

            r_hists_1.append(ROOT.TH2D('','', bin_settings[3][r_i], 0, 2*np.pi, len(bin_settings[2]), 0, bin_settings[0]))
            df_newer_1 = df_new_1[ df_new_1['r'] == r_i ]
            for i in range(len(df_newer_1)):
                r_hists_1[r_i].Fill(df_newer_1['tmid'].values[i], df_newer_1['rmid'].values[i],
                                df_newer_1['ly'].values[i] )

            r_hists_2.append(ROOT.TH2D('','', bin_settings[3][r_i], 0, 2*np.pi, len(bin_settings[2]), 0, bin_settings[0]))
            df_newer_2 = df_new_2[ df_new_2['r'] == r_i ]
            for i in range(len(df_newer_2)):
                r_hists_2[r_i].Fill(df_newer_2['tmid'].values[i], df_newer_2['rmid'].values[i],
                                df_newer_2['ly'].values[i] )

            r_hists_diff.append(ROOT.TH2D('','', bin_settings[3][r_i], 0, 2*np.pi, len(bin_settings[2]), 0, bin_settings[0]))
            df_newer_diff = df_new_diff[ df_new_diff['r'] == r_i ]
            for i in range(len(df_newer_diff)):
                r_hists_diff[r_i].Fill(df_newer_diff['tmid'].values[i], df_newer_diff['rmid'].values[i],
                                df_newer_diff['ly'].values[i] )


        z_hists_1.append(r_hists_1)
        z_hists_2.append(r_hists_2)
        z_hists_diff.append(r_hists_diff)

        if z_i < 3:
            padno = 2*z_i + 1
        elif z_i < 6:
            padno = 2*z_i + 7
        elif z_i < 9:
            padno = 2*z_i + 13
        elif z_i == 9:
            padno = 2*z_i + 19

        pad1 = c1.cd(padno)
        pad1.SetBorderSize(20)
        pad1.SetBorderMode(1)

        dummyH_list[0][z_i].Draw('colz')
        dummyH_list[0][z_i].SetTitle("Xerawdp") #%.2fcm < z < %.2fcm" %(z_i*3.03, (z_i+1)*3.03 
        dummyH_list[0][z_i].GetZaxis().SetTitle("<s1Area>")
        dummyH_list[0][z_i].GetXaxis().SetTitle("x position [cm]")
        dummyH_list[0][z_i].GetXaxis().SetTitleSize(0.06)
        dummyH_list[0][z_i].GetXaxis().CenterTitle()
        dummyH_list[0][z_i].GetYaxis().SetTitle("y position [cm]")
        dummyH_list[0][z_i].GetYaxis().SetTitleSize(0.06)
        dummyH_list[0][z_i].GetYaxis().CenterTitle()
        c1.SetRightMargin(0.2)

        for i in range(len(z_hists_1[z_i])):

            z_hists_1[z_i][i].GetZaxis().SetRangeUser(0, 7)
            z_hists_1[z_i][i].GetZaxis().SetTitle(peak + " ly [pe/keV]")
            z_hists_1[z_i][i].GetZaxis().SetTitleOffset(1.8)
            z_hists_1[z_i][i].Draw('pol colz a same')

        pad2 = c1.cd(padno+1)
        pad2.SetBorderSize(20)
        pad2.SetBorderMode(1)

        dummyH_list[1][z_i].Draw('colz')
        dummyH_list[1][z_i].SetTitle("Pax")
        dummyH_list[1][z_i].GetZaxis().SetTitle("<s1Area>")
        dummyH_list[1][z_i].GetXaxis().SetTitle("x position [cm]")
        dummyH_list[1][z_i].GetXaxis().SetTitleSize(0.06)
        dummyH_list[1][z_i].GetXaxis().CenterTitle()
        dummyH_list[1][z_i].GetYaxis().SetTitle("y position [cm]")
        dummyH_list[1][z_i].GetYaxis().SetTitleSize(0.06)
        dummyH_list[1][z_i].GetYaxis().CenterTitle()
        c1.SetRightMargin(0.2)

        for i in range(len(z_hists_2[z_i])):

            z_hists_2[z_i][i].GetZaxis().SetRangeUser(0, 7)
            z_hists_2[z_i][i].GetZaxis().SetTitle(peak + " ly [pe/keV]")
            z_hists_2[z_i][i].GetZaxis().SetTitleOffset(1.8)
            z_hists_2[z_i][i].Draw('pol colz a same')

        pad3 = c1.cd(padno+6)
        pad3.SetBorderSize(20)
        pad3.SetBorderMode(1)


        dummyH_list[2][z_i].Draw('colz')
        dummyH_list[2][z_i].SetTitle("|Difference|")
        dummyH_list[2][z_i].GetZaxis().SetTitle("<s1Area>")
        dummyH_list[2][z_i].GetXaxis().SetTitle("x position [cm]")
        dummyH_list[2][z_i].GetXaxis().SetTitleSize(0.06)
        dummyH_list[2][z_i].GetXaxis().CenterTitle()
        dummyH_list[2][z_i].GetYaxis().SetTitle("y position [cm]")
        dummyH_list[2][z_i].GetYaxis().SetTitleSize(0.06)
        dummyH_list[2][z_i].GetYaxis().CenterTitle()
        c1.SetRightMargin(0.2)

        for i in range(len(z_hists_diff[z_i])):

            z_hists_diff[z_i][i].GetZaxis().SetRangeUser(0, 1)
            z_hists_diff[z_i][i].GetZaxis().SetTitle(peak + " ly [pe/keV]")
            z_hists_diff[z_i][i].GetZaxis().SetTitleOffset(1.8)
            z_hists_diff[z_i][i].Draw('pol colz a same')
    
    for i in range(1,49):
        pad = c1.cd(i)
        pad.SetTopMargin(0.2)
        pad.SetBottomMargin(0)
    
    c1.Print(figure_name)
    c1.Clear()
    return
