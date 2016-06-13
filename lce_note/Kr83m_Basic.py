import hax

class Kr83m_Basic(hax.minitrees.TreeMaker):
    
    __version__ = '0.0.2'
    extra_branches = ['dataset_name','peaks.n_contributing_channels','peaks.hit_time_mean','peaks.left','peaks.right']
    
    def extract_data(self, event):
        
        # If there are no interactions at all, we can't extract anything...
        if not len(event.interactions):
            return dict()
        
        # Extract dataset_number
        dsetname = event.dataset_name
        if dsetname.endswith('.xed'):
            filename = dsetname.split("/")[-1]
            _, date, time, _ = filename.split('_')
            dataset_number = int(date) * 1e4 + int(time)
        else:
            # TODO: XENON1T support
            dataset_number = 0
            
        event_data = dict(event_number=event.event_number,
                          event_time=event.start_time,
                          dataset_number=dataset_number)
        
        # shortcuts for pax classes
        peaks = event.peaks
        interactions = event.interactions
        
        # assume 1st kr signal is interactions[0]
        s10 = interactions[0].s1
        s20 = interactions[0].s2
        
        # find 2nd kr interaction
        krInt = [0,0]  
        for i, interaction in enumerate(interactions):
            if interaction.s1 != s10 and interaction.s2 == s20 and krInt[0] == 0:
                krInt[0] = i
            elif interaction.s1 != s10 and interaction.s2 != s20 and krInt[1] == 0:
                krInt[1] = i
    
        # Distinction b/w single and double s2 events
        # Cut events without second s1
        if krInt[1] != 0:
            s11 = interactions[krInt[1]].s1
            s21 = interactions[krInt[1]].s2
            sInt = krInt[1]
        elif krInt[0] != 0:
            s11 = interactions[krInt[0]].s1
            s21 = -1
            sInt = krInt[0]
        else: return dict()
                
        # Find additional unwanted peaks
        s12 = -1
        for s1 in event.s1s:
            if s1 not in [s10,s11]: 
                s12 = s1
                break
        s22 = -1
        for s2 in event.s2s:
            if s2 not in [s20,s21]:
                s22 = s2
                break
                
        ##### Grab Data #####
    
        event_data.update(dict( s10Area = peaks[s10].area,
                                cs10Area = peaks[s10].area*interactions[0].s1_area_correction,
                                s10Coin = peaks[s10].n_contributing_channels,
                                s10Time = peaks[s10].hit_time_mean,
                                s10LeftEdge = peaks[s10].left*10.0,
                                s10RightEdge = peaks[s10].right*10.0,
                                s10x = interactions[0].x,
                                s10y = interactions[0].y,
                                s10z = interactions[0].z,
                                s20Area = peaks[s20].area,
                                cs20Area = peaks[s20].area*interactions[0].s2_area_correction,
                                s20Coin = peaks[s20].n_contributing_channels,
                                s20Time = peaks[s20].hit_time_mean,
                                s20LeftEdge = peaks[s20].left*10.0,
                                s20RightEdge = peaks[s20].right*10.0,
                                s11Area = peaks[s11].area,
                                cs11Area = peaks[s11].area*interactions[sInt].s1_area_correction,
                                s11Coin = peaks[s11].n_contributing_channels,
                                s11Time = peaks[s11].hit_time_mean,
                                s11LeftEdge = peaks[s11].left*10.0,
                                s11RightEdge = peaks[s11].right*10.0,
                                s11x = interactions[sInt].x,
                                s11y = interactions[sInt].y,
                                s11z = interactions[sInt].z ))
    
        if s21 != -1:
            s21Area = peaks[s21].area
            s21Coin = peaks[s21].n_contributing_channels
            s21Time = peaks[s21].hit_time_mean
            s21LeftEdge = peaks[s21].left*10.0
            s21RightEdge = peaks[s21].right*10.0
        else:
            s21Area = 0
            s21Coin = 0
            s21Time = 0
            s21LeftEdge = 0
            s21RightEdge = 0
    
        if s12 != -1:
            s12Area = peaks[s12].area
            s12Coin = peaks[s12].n_contributing_channels
            s12Time = peaks[s12].hit_time_mean
        else:
            s12Area = 0
            s12Coin = 0
            s12Time = 0
        
        if s22 != -1:
            s22Area = peaks[s22].area
            s22Coin = peaks[s22].n_contributing_channels
            s22Time = peaks[s22].hit_time_mean
        else:
            s22Area = 0
            s22Coin = 0
            s22Time = 0
            
        event_data.update(dict( s21Area = s21Area,
                                s21Coin = s21Coin,
                                s21Time = s21Time,
                                s21LeftEdge = s21LeftEdge,
                                s21RightEdge = s21RightEdge,
                                s12Area = s12Area,
                                s12Coin = s12Coin,
                                s12Time = s12Time,
                                s22Area = s22Area,
                                s22Coin = s22Coin,
                                s22Time = s22Time ))
        
        return event_data        
