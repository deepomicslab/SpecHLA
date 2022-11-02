import vcf
import graph

class Record:
    def __init__(self, pos = 0, gt0 = 0, gt1 = 0, PS = 0, idx=0):
        self.pos = pos
        if gt0 == 2 or gt1 == 2: # genotype 1/2: 1 -> 0 | 2 -> 1
            gt0 -= 1
            gt1 -= 1
        self.hap = gt0
        self.hap1 = gt1 
        self.ps = PS
        self.idx = idx

    def copy_from_rec(self, rec: vcf.model._Record, PS: int, idx:int):
        self.pos = rec.POS
        gt_str = rec.samples[0]['GT']
        gt0 = int(gt_str[0])
        gt1 = int(gt_str[2])
        if gt0 == 2 or gt1 == 2: # genotype 1/2: 1 -> 0 | 2 -> 1
            gt0 -= 1
            gt1 -= 1
        self.hap = gt0
        self.ps = PS
        self.idx = idx
    
        
    def flip(self):
        if self.hap == 0:
            self.hap = 1
        if self.hap == 1:
            self.hap = 0
    
    def phased(self):
        return self.ps != 0

    def copy_info(self, record):
        if self.pos != record.pos:
            return
        self.hap = record.pos
        self.ps = record.ps

    def switched(self, adjacent_record):
        if self.hap == adjacent_record.hap:
            return False
        else:
            return True
    
    def finalize_record(self, rec: vcf.model._Record):
        if not self.phased():
            return

        gt_str = rec.samples[0]['GT']
        gt0 = int(gt_str[0])
        gt1 = int(gt_str[2])
        hap0 = self.hap
        hap1 = abs(self.hap - 1)
        if gt0 == 2 or gt1 == 2:
            hap0 += 1
            hap1 += 1
        
        gt_str = str(hap0) + '|' + str(hap1)

        if self.phased():
            if ('PS' in rec.FORMAT.split(':')):
                rec.samples[0].data = rec.samples[0].data._replace(GT=gt_str, PS=self.ps)
                
            else:
                rec.add_format('PS')
                samp_fmt = vcf.model.make_calldata_tuple(rec.FORMAT.split(':'))
                tmp = rec.samples[0].data._asdict()
                tmp['PS'] = self.ps
                tmp['GT'] = gt_str
                rec.samples[0].data = samp_fmt(**tmp)

            rec.samples[0].gt_nums = gt_str
        


class PhaseSet:
    def __init__(self, starting_pos):
        self.starting_pos = starting_pos
        self.records_idx = set()
        self.records = dict()

    def contain_record(self, idx):
        return idx in self.records_idx
    
    def insert_record(self, record : Record):
        self.records_idx.add(record.pos)
        self.records[record.pos] = record

    def flip(self):
        for record in self.records.values():
            record.flip()
    
    def finalize_phaseset_label(self):
        self.starting_pos = min(self.records_idx) 
        for record in self.records.values():
            record.ps = self.starting_pos
    

    def intersect_phase_set(self, phase_set):
        intersection_length = 0
        hap0_support_length = 0
        extra_record_id = set()
        s_record: Record
        for s_record in phase_set.records:
            if s_record.pos in self.records_idx:
                intersection_length += 1
                record = self.records[s_record.pos]
                if record.hap == s_record.hap:
                    hap0_support_length += 1
            else:
                extra_record_id.add(s_record.pos)
        need_flip = False
        if intersection_length != 0:
            if hap0_support_length < intersection_length / 2:
                need_flip = True
        
        return need_flip, extra_record_id


class ChromosomoHaplotype:
    def __init__(self, in_vcf: vcf.Reader, chromo: str):
        self.chromo_record = dict()
        self.chromo_phase_set = dict()
        self.chromo_record2phaseset_map = dict()
        self.graph_struct = graph.Graph()
        rec:vcf.model._Record
        ps_label_fix = dict()
        idx = 0
        
        for rec in in_vcf.fetch(chromo):
            het = rec.samples[0].gt_type
            if het != 1:        # not het loci
                continue
            PS_fix = 0
            if rec.samples[0].phased:
                fmt = rec.FORMAT.split(':')
                if 'PS' in fmt:
                    PS = rec.samples[0]['PS']
                    if PS in ps_label_fix.keys():
                        PS_fix = ps_label_fix[PS]
                    else:
                        ps_label_fix[PS] = rec.POS
                        PS_fix = rec.POS
                else:
                    PS_fix = 1
            record = Record()
            record.copy_from_rec(rec, PS_fix, idx)
            idx += 1
            self.chromo_record[record.pos] = record
            if record.ps != 0:
                PS = record.ps
                self.chromo_record2phaseset_map[record.pos] = PS
                phase_set : PhaseSet
                if PS in self.chromo_phase_set.keys():
                    phase_set = self.chromo_phase_set[PS]
                else:
                    phase_set = PhaseSet(record.ps)
                    self.chromo_phase_set[PS] = phase_set
                phase_set.insert_record(record)            


    def construct_connection_graph(self, chromo_haplotype):
        """construct_connection_graph [create a graph structure based on primary and secondary phase set]
            node: phase set
            edge: shared record between phaseset 

        Args:
            chromo_haplotype (ChromosomoHaplotype): [description]
        """

        for phase_set in self.chromo_phase_set.values():
            self.graph_struct.insert_node(graph.Node(phase_set.starting_pos, 0)) # primary node
        for phase_set in chromo_haplotype.chromo_phase_set.values():
            self.graph_struct.insert_node(graph.Node(phase_set.starting_pos, 1)) # secondary Node

        phase_set : PhaseSet
        for phase_set in self.chromo_phase_set.values():
            met_phase_set = set()
            for record_pos in phase_set.records_idx:
                if record_pos not in chromo_haplotype.chromo_record2phaseset_map.keys(): # not in secondary phase set
                    continue
                ps_secondary = chromo_haplotype.chromo_record2phaseset_map[record_pos]
                if ps_secondary in met_phase_set: # connection already found
                    continue
                self.graph_struct.add_edge_with_ps(phase_set.starting_pos, ps_secondary)
            met_phase_set.clear()

        self.graph_struct.load_connected_components()     
        

    def merge_chromo_haplotype(self, chromo_haplotype):
        """merge_chromo_haplotype [merge chromosome level haplotype based on the created graph structure]
            note that connected component are listed in such order: s - f -s or f-s-f
        Args:
            chromo_haplotype (ChromosomoHaplotype): [secondary chromosome level hap]
        """
        for connected_phase_sets in self.graph_struct.connected_components:
            n_nodes = len(connected_phase_sets)
            if n_nodes == 1:
                node = self.graph_struct.get_node(connected_phase_sets[0])
                if node.is_primary():   # no new info
                    continue
                secondary_phase_set = chromo_haplotype.chromo_phase_set[node.ps]
                primary_phase_set = PhaseSet(secondary_phase_set.starting_pos)                               
                self.create_phase_set_from_secondary(primary_phase_set, secondary_phase_set)
            
            else:   # always start from a primary phase set, then BFS
                start_node_id = self.graph_struct.get_closeset_primary_node(connected_phase_sets[0])
                start_phase_set = self.chromo_phase_set[self.graph_struct.get_node(start_node_id).ps]
                visited = dict()
                for id in connected_phase_sets:
                    visited[id] = False
                queue = []
                queue.append(start_node_id)
                visited[start_node_id] = True

                while queue:
                    s = queue.pop(0)
                    if s == start_node_id:
                        continue
                    node = self.graph_struct.get_node(s)
                    if node.is_primary():
                        r_phase_set = self.chromo_phase_set[node.ps]
                        self.connect_phase_set(start_phase_set, r_phase_set)
                    else: 
                        secondary_phase_set = chromo_haplotype.chromo_phase_set[node.ps]
                        self.extend_phase_set(start_phase_set, secondary_phase_set)
                    for i in self.graph_struct.adj_list[s]: 
                        if visited[i] == False: 
                            queue.append(i) 
                            visited[i] = True
    
    
    def create_phase_set_from_secondary(self, phase_set: PhaseSet, secondary_phase_set: PhaseSet):
        """create_phase_set_from_secondary [new phase set introduced from secondary, no intersection with existing phase set]

        Args:
            phase_set (PhaseSet): [primary phase set]
            secondary_phase_set (PhaseSet): [secondary phase set]
        """
        for secondary_record in secondary_phase_set.records.values():
            # this means potential bug, should not be here
            if secondary_record.pos in self.chromo_record2phaseset_map.keys():
                continue
            record = self.chromo_record[secondary_record.pos]
            record.copy_info(secondary_record)
            phase_set.insert_record(record)
            self.chromo_record2phaseset_map[record.pos] = phase_set.starting_pos
        self.chromo_phase_set[phase_set.starting_pos] = phase_set


    def extend_phase_set(self, phase_set: PhaseSet, secondary_phase_set: PhaseSet):
        """extend_phase_set [merge two phaseset with intersection or superset relationship]

        Args:
            phase_set (PhaseSet): [description]
            secondary_phase_set (PhaseSet): [description]
        """
        need_flip, extra_record_pos = phase_set.intersect_phase_set(secondary_phase_set)
        if len(extra_record_pos) == 0: # super set condition, do nothing
            return
        if need_flip:
            secondary_phase_set.flip()
        for record_pos in extra_record_pos:
            if record_pos not in self.chromo_record2phaseset_map.keys(): # simple extend the set
                record = self.chromo_record[record_pos]
                record.copy_info(secondary_phase_set.records[record_pos])
                self.chromo_record2phaseset_map[record_pos] = phase_set.starting_pos
                phase_set.insert_record(record)
            else:
                record = secondary_phase_set.records[record_pos]
                phase_set.insert_record(record)
            
    
    
    def connect_phase_set(self, f_phase_set: PhaseSet, s_phase_set:PhaseSet):
        need_flip, extra_record_pos = f_phase_set.intersect_phase_set(s_phase_set)
        if need_flip:
            s_phase_set.flip()
        for record_pos, record in s_phase_set.records:
            f_phase_set.insert_record(record)
            self.chromo_record2phaseset_map[record_pos] = f_phase_set.starting_pos
        
        del self.chromo_phase_set[s_phase_set.starting_pos]
        s_phase_set.records.clear()
    


    def finalize_chromosome_haplotype(self):
        self.graph_struct.clear()
        self.chromo_record2phaseset_map.clear()
        phase_set: PhaseSet
        temp = dict()
        for phase_set in self.chromo_phase_set.values():
            phase_set.finalize_phaseset_label()
            temp[phase_set.starting_pos] = phase_set
            for record_idx in phase_set.records_idx:
                self.chromo_record2phaseset_map[record_idx] = phase_set.starting_pos
        


    