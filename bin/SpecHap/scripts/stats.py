class PhaseSetStats:
    def __init__(self, switch_count, mismatch_count, phased_snp, N50, AN50, spanned_snp):
        self.switch_count = switch_count
        self.mismatch_count = mismatch_count
        self.phased_snp = phased_snp
        self.N50 = N50
        self.AN50 = AN50
        self.spanned_snp = spanned_snp
    
    def get_AN50(self):
        return self.AN50
    
    def get_N50(self):
        return self.N50

    def get_phased_snp(self):
        return self.phased_snp
    
    def get_poss_mismatch(self):
        return self.phased_snp

    def get_poss_switch(self):
        if self.phased_snp <= 3:
            return 0
        else:
            return self.phased_snp - 3
    
    def get_mismatch_count(self):
        return self.mismatch_count
    
    def get_switch_count(self):
        if self.switch_count < 0:
            return 0
        return self.switch_count

    def get_switch_error(self):
        poss_sw = self.get_poss_switch()
        if poss_sw < 1:
            return 0
        else:
            return self.get_switch_count()  / poss_sw
    
    def get_mismatch_error(self):
        poss_mm = self.get_poss_mismatch()
        if poss_mm < 1 :
            return 0
        else:
            return self.get_mismatch_count() / poss_mm
    
    def get_spanned_snp(self):
        return self.spanned_snp
    
    

class HapStats:
    def __init__(self, total_snp_count = 0, total_base_pair_count = 0):
        self.phase_set_stats = list()
        self.total_snp_count = total_snp_count
        self.total_base_pair_count = total_base_pair_count

    def insert_phase_set_stats(self, idx:int, phase_set_stats: PhaseSetStats):
        self.phase_set_stats.append(phase_set_stats)
    
    def insert_hap_stats(self, hapstats):
        self.total_snp_count += hapstats.total_snp_count
        self.total_base_pair_count += hapstats.total_base_pair_count
        for i in hapstats.phase_set_stats:
            self.phase_set_stats.append(i)

    def get_total_phased(self):
        total = 0
        for phase_set_stats in self.phase_set_stats:
            if phase_set_stats.get_phased_snp() < 2: 
                continue
            total += phase_set_stats.get_phased_snp()
        return total

    def get_AN50(self):
        AN50 = 0
        phase_sets = self.phase_set_stats
        phase_sets.sort(key=lambda x: x.AN50, reverse=True)
        phased_sum = 0
        for phase_set_stats in phase_sets:
            phased_sum += phase_set_stats.get_phased_snp()
            if phased_sum > self.total_snp_count / 2:
                AN50 = phase_set_stats.get_AN50()
                break
        return AN50

    def get_N50(self):
        N50 = 0
        phase_sets = self.phase_set_stats
        phase_sets.sort(key=lambda x: x.N50, reverse=True)
        total_span = 0
        for phase_set_stats in phase_sets:
            total_span += phase_set_stats.get_N50()
            if total_span > self.total_base_pair_count / 2:
                N50 = phase_set_stats.get_N50()
                break
        return N50
    
    def get_mismatch_error(self):
        mismatch_count = 0
        poss_mismatch = 0
        phase_set_stats : PhaseSetStats
        for phase_set_stats in self.phase_set_stats:
            if phase_set_stats.get_poss_mismatch() < 1:
                continue
            mismatch_count += phase_set_stats.get_mismatch_count()
            poss_mismatch += phase_set_stats.get_poss_mismatch()
        
        if poss_mismatch == 0:
            return 0
        return mismatch_count / poss_mismatch
    
    def get_switch_error(self):
        switch_count = 0
        poss_switch = 0
        phase_set_stats : PhaseSetStats
        for phase_set_stats in self.phase_set_stats:
            if phase_set_stats.get_poss_switch() < 1:
                continue
            switch_count += phase_set_stats.get_switch_count()
            #if phase_set_stats.get_switch_count() != 0:
                #print (phase_set_stats.AN50)
            poss_switch += phase_set_stats.get_poss_switch()
        
        if poss_switch == 0:
            return 0
        
        #print(switch_count)
        return switch_count / poss_switch
    
    
    def get_total_spanned(self):
        total = 0
        for phase_set_stats in self.phase_set_stats:
            total += phase_set_stats.spanned_snp
        
        return total
