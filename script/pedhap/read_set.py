from pedigree import Trio


class Read:
    def __init__(
            self,
            mapq: int,
            block_id: int,
            threshold1: float = 0.1,
            threshold2: float = 0,
    ):
        self.block_id = block_id
        self.mapq = mapq
        self.covered_blocks = {}
        self.blocks = []
        self.block_reverses = []
        self.uncertain_blocks = []
        self.support_situation = {}
        self.confilict_side = {}
        self.threshold1 = threshold1
        self.threshold2 = threshold2

    def get_confilict_poses(self):
        res = []
        for k in self.confilict_side.keys():
            res = res + self.support_situation[k][self.confilict_side[k]]
        return res

    def set_covered_block(self, b_id: int, side: int, pos: int, value=1):
        # self.support_situation[side].append(pos)
        if b_id in self.covered_blocks.keys():
            item = self.covered_blocks[b_id]
        else:
            self.support_situation[b_id] = {0:[],1:[]}
            # self.confilict_side[b_id] = {}
            self.covered_blocks[b_id] = [0, 0]
            item = self.covered_blocks[b_id]
        self.support_situation[b_id][side].append(pos)
        if side == 0:
            item[0] = item[0] + value
        else:
            item[1] = item[1] + value

    def init_blocks(self, ensure_block = None):
        # if len(self.covered_blocks) <2:
        #     return
        for k, v in sorted(self.covered_blocks.items()):
            if ensure_block != None and len(ensure_block) != 0 and k == ensure_block[0]:
                if ensure_block[1] == 0:
                    self.confilict_side[k] = 1
                    self.block_reverses.append(False)
                    self.blocks.append(k)
                    continue
            if v[1] == v[0] or (v[1] != 0 and v[0] != 0 and abs(max(v)/min(v)) <= self.threshold1) or ( abs(v[1]- v[0]) <= self.threshold2 and (v[1] == 0 or v[0] == 0)):
                self.uncertain_blocks.append(k)
            else:
                need_reverse = False
                self.confilict_side[k] = 1
                if v[1] > v[0]:
                    need_reverse = True
                    self.confilict_side[k] = 0
                self.blocks.append(k)
                self.block_reverses.append(need_reverse)
        self.blocks.append(self.block_id)
        self.block_reverses.append(False)

    def get_blocks_info(self):
        return self.blocks, self.block_reverses, self.uncertain_blocks


class ReadSet(object):
    def __init__(self):
        self.reverse_info = {}
        self.father_dict = {}
        self.size_dict = {}
        self.uncertain_blocks = []
        self.confilict_poses = []

    def contains_phasing_info(self):
        tag = False
        for k, v in self.size_dict.items():
            if v > 2:
                tag = True
        return tag

    def get_phase_id(self, block_id):
        if block_id not in self.reverse_info.keys():
            return (0, 0)
        else:
            return (self.find(block_id), self.reverse_info[block_id])

    def add_read(self, read: Read, ensure_block = None):
        read.init_blocks(ensure_block)
        block_ids, reverses, uncertain_blocks = read.get_blocks_info()
        print(block_ids, "ccss")
        self.confilict_poses = self.confilict_poses + read.get_confilict_poses()
        for b in uncertain_blocks:
            if b not in self.uncertain_blocks:
                self.uncertain_blocks.append(b)
        # 如果第一个block已经存在并且有冲突
        first_block = block_ids[0]
        if first_block in self.reverse_info.keys():
            if reverses[0] != self.reverse_info[first_block]:
                for i in range(0, len(reverses)):
                    reverses[i] = not reverses[i]
        # 记录block是否反转并且加入uf
        for i in range(0, len(block_ids)):
            b_id = block_ids[i]
            print(block_ids,"3333")
            r = reverses[i]
            if b_id not in self.father_dict.keys():
                self.father_dict[b_id] = b_id
                self.size_dict[b_id] = 1
            self.reverse_info[b_id] = r
            self.union(b_id, first_block)

    # def add_block(self, block_id: int, need_reverse: bool):
    #     node = Node(block_id, need_reverse)
    #     if node in self.nodes_set:

    def find(self, node):
        father = self.father_dict[node]
        if (node != father):
            if father != self.father_dict[father]:  # 在降低树高优化时，确保父节点大小字典正确
                self.size_dict[father] -= 1
            father = self.find(father)
        self.father_dict[node] = father
        return father

    def is_same_set(self, node_a, node_b):
        return self.find(node_a) == self.find(node_b)

    def union(self, node_a, node_b):
        if node_a == node_b:
            return
        if node_a is None or node_b is None:
            return
        a_head = self.find(node_a)
        b_head = self.find(node_b)

        if (a_head != b_head):
            a_set_size = self.size_dict[a_head]
            b_set_size = self.size_dict[b_head]
            if (a_set_size >= b_set_size):
                self.father_dict[b_head] = a_head
                self.size_dict[a_head] = a_set_size + b_set_size
            else:
                self.father_dict[a_head] = b_head
                self.size_dict[b_head] = a_set_size + b_set_size
