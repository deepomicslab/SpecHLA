class Node:
    def __init__(self, ps: int, categ: int):
        self.ps = ps
        self.categ = categ
        self.id = -1
        return
    
    def is_primary(self):
        return self.categ == 0

    def is_secondary(self):
        return self.categ == 1
    
    def set_id(self, id: int):
        self.id = id


class Graph:
    def __init__(self):
        self.n_node = 0
        self.nodes = dict()
        self.adj_list = dict()
        self.connected_components = list()
        self.primary_node_ps_id_map = dict()
        self.secondary_node_ps_id_map = dict()

    def clear(self):
        self.n_node = 0
        self.nodes.clear()
        self.adj_list.clear()
        self.connected_components.clear()
        self.primary_node_ps_id_map.clear()
        self.secondary_node_ps_id_map.clear()
        
    def insert_node(self, node : Node):
        node.set_id(self.n_node)
        self.nodes[self.n_node] = node
        if node.is_primary():
            self.primary_node_ps_id_map[node.ps] = self.n_node
        if node.is_secondary():
            self.secondary_node_ps_id_map[node.ps] = self.n_node 
        self.n_node += 1
    
    def add_edge(self, s:int, e:int):
        if s not in self.adj_list.keys():
            self.adj_list[s] = list()
        self.adj_list[s].append(e)

        if e not in self.adj_list.keys():
            self.adj_list[e] = list()
        self.adj_list[e].append(s)

    def get_node_id(self, ps: int, is_primary:bool):
        if is_primary:
            return self.primary_node_ps_id_map[ps]
        else:
            return self.secondary_node_ps_id_map[ps]

    def get_node(self, id: int):
        return self.nodes[id]
            

    def get_closeset_primary_node(self, id:int):
        node = self.get_node(id)
        if node.is_primary():
            return id
        else: 
            return self.adj_list[id][0]

    def add_edge_with_ps(self, s_ps, e_ps):
        s = self.get_node_id(s_ps, True)
        e = self.get_node_id(e_ps, False)
        self.add_edge(s, e)
    
    def DFSUtil(self, temp:list, node_id:int, visited:list): 
        visited[node_id] = True
        temp.append(node_id)
        for i in self.adj_list[node_id]: 
            if visited[i] == False: 
                temp = self.DFSUtil(temp, i, visited) 
        return temp 

    def load_connected_components(self):
        if len(self.connected_components) == 0:
            self.connected_components.clear()
        visited = [False for i in range(0, self.n_node)]
        for i in range(0, self.n_node):
            if i not in self.adj_list.keys():
                self.adj_list[i] = list()
        for node_id in range(self.n_node): 
            if visited[node_id] == False: 
                temp = [] 
                self.connected_components.append(self.DFSUtil(temp, node_id, visited))
        