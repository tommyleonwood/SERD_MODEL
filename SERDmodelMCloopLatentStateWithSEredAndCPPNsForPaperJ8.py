import numpy as np 
import random
import math
import copy


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from tqdm import trange, tqdm
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from sklearn.cluster import KMeans
import networkx as nx
from collections import deque

from matplotlib.animation import FFMpegWriter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec
from matplotlib.animation import FuncAnimation
from scipy.spatial import cKDTree
from matplotlib.widgets import Slider


# ────────────────────────────────────────────────────────────────────
#  SERD: user–visible run-time parameters
#  (edit these; everything further down just reads them)
# ────────────────────────────────────────────────────────────────────

# ---------------- Ensemble size & run length ------------------------
M          = 10   # number of independent Monte-Carlo replicas
T_max      = 30     # SERD time-steps per replica

# ---------------- Warm-up / initial condition -----------------------
use_latent     = False   # True → latent warm-up, False → straight chain
latent_steps   = 20     # ticks for the warm-up run if use_latent=True
init_SE_filament = 10   # straight-chain length if use_latent=False

# ---------------- Stochastic rule probabilities ---------------------
pp_split_probability = 0.50  # P(PP → two PPs) at t = 0
duplication_probability = 1/3  # P(SE duplication) per idle SE
reduction_probability   = 1/3  # P(SE reduction)  per idle SE
decay_gamma   = 1    # exponential decay of pp_split_probability
t_stop_splits = 10       # force p_split = 0 beyond this global tick

# ---------------- Embedding / visualisation -------------------------
embed_freq   = 0.1  # chance of doing one spring-electric step each tick
axisScale    = init_SE_filament  # axis limits in 3-D scatter plot

# Hooke-spring constants & Coulomb repulsion
k_springSE       = 1.0   # IG–IG spring (along SE)
k_springPPtoIG   = 2.0   # PP tether to its IG neighbours
k_repulsionIn    = 1.5   # IG–IG repulsion strength
repulsion_radiusIn = 10  # cut-off radius for repulsion

# Misc. embedding tweaks
SE_lengthIn        = 1.0   # natural length of IG–IG spring
PPtoIGlength       = 0.1   # natural length of PP–IG spring
randomOffsetIn     = 0.05  # thermal jiggle added after each embed step
maxRandomOffsetIn  = 0.10  # hard cap on that jiggle
minNodeUpdateIn    = 0.01  # ignore displacements smaller than this

# ---------------- Display styling (size/opacity only) ---------------
pp_alpha   = min(1.0, 100.0/M)  # PP dot opacity
ig_alpha   = min(1.0,  10.0/M)  # IG dot opacity
edge_alpha = min(0.5,   2.0/M)  # line opacity
pp_size    = 10                 # PP dot size (points)
ig_size    = 5                  # IG dot size
edge_width = 0.001              # line width in scatter plot

# ---------------- Static colours & symbols --------------------------
initial_edge_color     = "blue"
duplicating_edge_color = "lightgreen"
remaining_edge_color   = "yellow"
updated_edge_color     = "darkgreen"
propagating_PSM_color  = "orange"
PP_node_color          = "black"
IG_node_color          = "lightgrey"
PP_edge_color          = "black"
fused_PP_edge_color    = "grey"
background_color       = "white"

# ---------------- Book-keeping globals (do not edit) ----------------
psm_queue, next_psm_queue = [], []  # double-buffer for PSM sweeps
next_psm_id  = 0                    # monotonically increasing PSM label
force_quivers = []                  # matplotlib quiver handles (debug)
nodeOffset    = 0.6                 # radial offset when spawning new IGs
simulation_halted = False           # emergency stop flag

# ---------------- Data structures (filled later) --------------------
IG_dict, SE_dict, PP_dict = {}, {}, {}
pp_ig_edges = {}            # PP → set(IG) incidence map

# ---------------- ID counters for fresh nodes/edges -----------------
next_pp_id, next_se_id, next_ig_id = 2, 1, 2

# ---------------- Histories recorded by the analysis code ----------
pp_history, ig_history = [], []   # will hold node positions per TS



# Initialize IG, SE, and PP dictionaries with explicit neighborhood definitions
IG_dict = {}

SE_dict = {}

PP_dict = {}

# The passive edge mapping between PPs and IGs:
pp_ig_edges = {}


next_pp_id=2
next_se_id=1
next_ig_id=2

# ─── NEW GLOBALS: decay schedule & history ──────────────────────────────────

# (a) total TS over which p_split goes from p_init → 0
inflation_duration = 5 
p_split_init       = 1.0  

# (b) histories of node‐positions per TS
pp_history = []   # will hold np.array of shape (N_pp,3) each TS
ig_history = []   # will hold np.array of shape (N_ig,3) each TS

# (c) helper to build a straight‐chain of length n between two PPs
def initialize_chain(n):
    """n SEs in a line, (n+1) IGs between two PPs at the ends."""
    global IG_dict, SE_dict, PP_dict, pp_ig_edges,next_se_id,next_ig_id
    IG_dict.clear(); SE_dict.clear(); PP_dict.clear(); pp_ig_edges.clear()

    backShift=math.ceil(n/2)
    # 1) make IGs at x=0..n
    for i in range(n+1):
        ig = f"n{i}"
        IG_dict[ig] = {"position": [float(i)-backShift,0,0], "neighborhood": set(), "psm_ids": set()}

    # 2) make SEs s0..s{n-1} linking n{i}→n{i+1}
    for i in range(n):
        se = f"s{i}"
        SE_dict[se] = {"nodes": (f"n{i}",f"n{i+1}")}
        IG_dict[f"n{i}"]["neighborhood"].add(se)
        IG_dict[f"n{i+1}"]["neighborhood"].add(se)

    # 3) place PPs outside, at x=-1 and x=n+1
    PP_dict["p0"] = {"position":[-1.0-backShift,0,0], "neighborhood": {"n0"}}
    PP_dict["p1"] = {"position":[n+1.0-backShift,0,0], "neighborhood": {f"n{n}"}}
    IG_dict["n0"]["neighborhood"].add("p0")
    IG_dict[f"n{n}"]["neighborhood"].add("p1")
    pp_ig_edges["p0"] = {"n0"}
    pp_ig_edges["p1"] = {f"n{n}"}

    next_se_id=n
    next_ig_id=n+1


print(PP_dict)
print(IG_dict)
print(SE_dict)


PSM_dict = {}
duplication_queue = list(SE_dict.keys())
updated_SEs = set()
new_SEs_after_split = []



def compute_offset(pos1: np.ndarray, pos2: np.ndarray, scale: float = 0.2) -> np.ndarray:
    direction = pos2 - pos1
    norm = np.linalg.norm(direction)
    if norm == 0:
        return np.array([0.1, 0.1, 0]) * scale
    return (direction / norm) * scale


def get_new_psm_id():
    global next_psm_id
    new_id = f"psm{next_psm_id}"
    next_psm_id += 1
    return new_id

def get_new_pp_id():
    global next_pp_id
    new_id = f"p{next_pp_id}"
    next_pp_id += 1
    return new_id

def get_new_se_id():
    global next_se_id
    new_id = f"s{next_se_id}"
    next_se_id += 1
    return new_id

def get_new_ig_id():
    global next_ig_id
    new_id = f"n{next_ig_id}"
    next_ig_id += 1
    return new_id




def draw_forces(pp_positions, masses):
    global force_quivers
    # remove old arrows
    for q in force_quivers:
        q.remove()
    force_quivers = []

    if len(pp_positions) < 2:          # nothing to do
        return

    # crude O(N^2) force calc – fine for ≤ 200 PPs
    for i, pi in enumerate(pp_positions):
        Fi = np.zeros(3)
        for j, pj in enumerate(pp_positions):
            if i == j: 
                continue
            r = pj - pi
            d = np.linalg.norm(r)
            if d == 0: 
                continue
            Fi += masses[i]*masses[j] * r / d**3
        # scale the arrow for plotting aesthetics
        arrow = ax.quiver(pi[0], pi[1], pi[2],
                          Fi[0], Fi[1], Fi[2],
                          length=0.5, normalize=True, color='crimson')
        force_quivers.append(arrow)


# --- Invariant Verification Function ---
def verify_psm_invariant(IG_dict, PSM_dict):
    """
    For every IG with exactly two PSMs, verify that the union of the FROM sets
    of those PSMs equals the IG’s neighborhood (or that the IG’s degree is equal to
    len(union(FROM)) + 1). If not, print an error.
    """
    invariant_ok = True
    for ig_id, ig in IG_dict.items():
        psm_ids = ig.get("psm_ids", set())
        if len(psm_ids) == 2:
            union_from = set()
            for psm_id in psm_ids:
                union_from |= PSM_dict.get(psm_id, {}).get("FROM", set())
            expected_degree = len(union_from) + 1
            actual_degree = len(ig.get("neighborhood", set()))
            if actual_degree != expected_degree:
                invariant_ok = False
    return invariant_ok


def test_se_reduction_invariant(final_ig, IG_dict, PSM_dict):
    """
    Test that for the final IG (produced by an SE reduction) that has exactly two PSMs,
    the union of the FROM sets of those PSMs equals the IG's neighborhood.
    
    If the invariant is violated, print detailed information for debugging.
    """
    if final_ig not in IG_dict:
        print(f"[TEST] ERROR: Final IG {final_ig} is not present in IG_dict!")
        return False

    ig = IG_dict[final_ig]
    psm_ids = ig.get("psm_ids", set())

    # Only run the invariant check if there are exactly two PSMs.
    if len(psm_ids) == 2:
        union_from = set()
        for pid in psm_ids:
            psm = PSM_dict.get(pid)
            if psm is None:
                continue
            union_from |= psm.get("FROM", set())
        neighborhood = ig.get("neighborhood", set())

        if union_from == neighborhood:
            return True
        else:
            for pid in psm_ids:
                psm = PSM_dict.get(pid, {})
            return False
    else:
        return True

def bypass_reduction_due_to_shared_edge(se_id, ig1, ig2, IG_dict, PSM_dict, SE_dict):
    """
    Returns True if both bounding IGs (ig1 and ig2) have at least one PSM with se_id in its FROM set.
    In that case, we bypass the reduction (i.e. do not merge the IGs) and mark the SE as darkgreen,
    and then remove its "pending_op" flag so that the update-phase sees it as updated.
    Otherwise, returns False so that reduction can proceed.
    """
    # Get the sets of PSM IDs for each IG.
    psm_ids_ig1 = IG_dict[ig1].get("psm_ids", set())
    psm_ids_ig2 = IG_dict[ig2].get("psm_ids", set())
    
    # Check if at least one PSM on each IG has se_id in its FROM set.
    has_from_on_ig1 = any(se_id in PSM_dict[psm]["FROM"] for psm in psm_ids_ig1)
    has_from_on_ig2 = any(se_id in PSM_dict[psm]["FROM"] for psm in psm_ids_ig2)
    
    if has_from_on_ig1 and has_from_on_ig2:
        if se_id in SE_dict:
            # NEW: Remove the pending_op flag so that this SE is no longer considered pending.
            if "pending_op" in SE_dict[se_id]:
                del SE_dict[se_id]["pending_op"]
        return True
    else:
        return False



def estimate_topological_dimension(IG_dict, SE_dict, n_roots=5, min_dist=1, max_dist=None):
    """
    Estimate the graph‐distance (fractal) dimension D of the IG network:
      IG_dict: mapping ig_id -> {...}
      SE_dict: mapping se_id -> {'nodes': (ig1, ig2), ...}
      n_roots: number of random starting IGs to average over
      min_dist: smallest ℓ to include in the fit
      max_dist: largest ℓ to include (defaults to diameter)
    Returns:
      D_mean : mean estimated dimension
      D_std  : standard deviation over roots
      per_root : list of individual D estimates
    """
    # 1) build adjacency
    adj = {ig: set() for ig in IG_dict}
    for se in SE_dict.values():
        u,v = se['nodes']
        if u in adj and v in adj:
            adj[u].add(v)
            adj[v].add(u)

    # choose roots
    all_igs = list(adj)
    if n_roots >= len(all_igs):
        roots = all_igs
    else:
        roots = random.sample(all_igs, n_roots)

    Ds = []
    for root in roots:
        # 2) BFS for graph‐distance
        dist = {root:0}
        q = deque([root])
        while q:
            u = q.popleft()
            for w in adj[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1
                    q.append(w)

        # 3) form cumulative counts N(ℓ)
        max_d = max(dist.values())
        counts = np.zeros(max_d+1, int)
        for d in dist.values():
            counts[d] += 1
        cum = np.cumsum(counts)  # cum[ℓ] = #nodes ≤ ℓ

        # decide fitting range
        ℓ = np.arange(len(cum))
        if max_dist is None:
            fit_max = max_d
        else:
            fit_max = min(max_d, max_dist)
        mask = (ℓ>=min_dist)&(ℓ<=fit_max)&(cum>1)
        logℓ = np.log(ℓ[mask])
        logN = np.log(cum[mask])
        if len(logℓ) < 2:
            continue
        # 4) linear fit: log N = D·log ℓ + c
        D, _ = np.polyfit(logℓ, logN, 1)
        Ds.append(D)

    D_mean = float(np.mean(Ds)) if Ds else 0.0
    D_std  = float(np.std(Ds))  if Ds else 0.0
    return D_mean, D_std, Ds


# Spring-electric embedding function -------------------------------------------------------------
def spring_electric_embedding_local(IG_dict, SE_dict, PP_dict, 
                                      damping=0.1,
                                      max_disp=0.5,
                                      min_disp=0.01,
                                      k_spring1=k_springSE,
                                      k_spring2=k_springPPtoIG,
                                      rest_length=SE_lengthIn,
                                      rest_lengthPPtoIG=PPtoIGlength,
                                      k_repulsion=k_repulsionIn,
                                      repulsion_radius=repulsion_radiusIn,
                                      randomOffset=randomOffsetIn,
                                      maxRandomOffset=maxRandomOffsetIn):
    """
    Compute the new positions of IGs and PPs by applying spring‐electric forces.
    Uses a damping factor and a force cap (max_disp) to prevent runaway motion.
    If the computed displacement is below min_disp, the node’s position is left unchanged.
    This version caches the IG and PP keys and builds a mapping for IG lookups.
    """
    # Cache IG keys and build an index lookup.
    ig_keys = list(IG_dict.keys())
    ig_index = {k: i for i, k in enumerate(ig_keys)}
    num_IGs = len(ig_keys)
    
    # Build PP keys and their positions.
    pp_keys = list(PP_dict.keys())
    
    # Build positions arrays for IGs and PPs.
    ig_positions = np.array([IG_dict[k]["position"] for k in ig_keys], dtype=np.float64)
    pp_positions = np.array([PP_dict[k]["position"] for k in pp_keys], dtype=np.float64)
    positions = np.concatenate((ig_positions, pp_positions), axis=0)
    
    forces = np.zeros_like(positions, dtype=np.float64)
    
    # Process spring forces for SEs.
    for se_id, se in SE_dict.items():
        node1, node2 = se["nodes"]
        try:
            idx1 = ig_index[node1]
            idx2 = ig_index[node2]
        except KeyError as e:
            raise KeyError(f"[ERROR] spring_electric_embedding_local: Could not find IG in SE {se_id}: {e}")
        p1, p2 = positions[idx1], positions[idx2]
        distance = np.linalg.norm(p2 - p1)
        direction = (p2 - p1) / distance if distance != 0 else 0
        force = k_spring1 * (distance - rest_length)
        forces[idx1] += force * direction
        forces[idx2] -= force * direction

    # Process spring forces for PPs connected to IGs.
    for i, pp_id in enumerate(pp_keys):
        # For PP nodes, their index starts after the IG nodes.
        idx_pp = num_IGs + i
        p1 = positions[idx_pp]
        for ig_id in PP_dict[pp_id]["neighborhood"]:
            if ig_id not in ig_index:
                raise KeyError(f"[ERROR] PP {pp_id} references missing IG {ig_id} in its neighborhood.")
            idx_ig = ig_index[ig_id]
            p2 = positions[idx_ig]
            distance = np.linalg.norm(p2 - p1)
            if distance > 0:
                direction = (p2 - p1) / distance
                force = k_spring2 * (distance - rest_lengthPPtoIG)
                forces[idx_pp] += force * direction
                forces[idx_ig] -= force * direction

    # Electric repulsion between IG nodes.
    tree = cKDTree(positions[:num_IGs])
    for i, pos in enumerate(positions[:num_IGs]):
        neighbors = tree.query_ball_point(pos, r=repulsion_radius)
        for j in neighbors:
            if i != j:
                p2 = positions[j]
                distance = np.linalg.norm(p2 - pos)
                if distance > 0:
                    direction = (p2 - pos) / distance
                    repulsion = k_repulsion / (distance ** 2)
                    forces[i] -= repulsion * direction

    # Update positions for IG nodes.
    for idx, node_id in enumerate(ig_keys):
        raw_disp = forces[idx] * damping
        disp_norm = np.linalg.norm(raw_disp)
        if disp_norm < min_disp:
            continue  # Skip update if change is too small.
        disp = raw_disp if disp_norm <= max_disp else raw_disp * (max_disp / disp_norm)
        new_pos = positions[idx] + disp
        # Add a small random perturbation scaled by the displacement.
        rand_off = disp_norm * 0.2
        rand_off = min(rand_off, maxRandomOffset)
        new_pos += np.array([random.uniform(-rand_off, rand_off) for _ in range(3)])
        IG_dict[node_id]["position"] = new_pos.tolist()

    # Update positions for PP nodes.
    for i, pp_id in enumerate(pp_keys):
        idx_pp = num_IGs + i
        raw_disp = forces[idx_pp] * damping
        disp_norm = np.linalg.norm(raw_disp)
        if disp_norm < min_disp:
            continue
        disp = raw_disp if disp_norm <= max_disp else raw_disp * (max_disp / disp_norm)
        new_pos = positions[idx_pp] + disp
        rand_off = disp_norm * 0.2
        rand_off = min(rand_off, maxRandomOffset)
        new_pos += np.array([random.uniform(-rand_off, rand_off) for _ in range(3)])
        PP_dict[pp_id]["position"] = new_pos.tolist()
    return



# Define UPDATE OPERATIONS:----------------------
def update_neighbors(element_id, remove_items, add_items, element_dict):
    # Instead of silently doing nothing if the element is missing, raise an error.
    if element_id not in element_dict:
        raise KeyError(f"update_neighbors: Element {element_id} not found in dictionary.")
    element = element_dict[element_id]
    if "nodes" in element:
        nodes = list(element["nodes"])
        # We assume the lengths of remove_items and add_items are matched.
        remove_list = list(remove_items)
        add_list = list(add_items)
        for i in range(len(nodes)):
            for j in range(len(remove_list)):
                if nodes[i] == remove_list[j]:
                    old_val = nodes[i]
                    nodes[i] = add_list[j]
        element["nodes"] = tuple(nodes)
    elif "neighborhood" in element:
        for item in remove_items:
            if item not in element["neighborhood"]:
                raise KeyError(f"update_neighbors: Cannot discard {item} from neighborhood of {element_id} because it is not present.")
            element["neighborhood"].discard(item)
        if add_items:
            element["neighborhood"].update(add_items)

def update_pp_connection(pp_id, old_igs, new_igs, PP_dict, pp_ig_edges):
    # For each IG that is being replaced, if it isn’t present then raise an error.
    for old_ig in old_igs:
        if old_ig not in PP_dict[pp_id]["neighborhood"]:
            raise KeyError(f"update_pp_connection: PP {pp_id} does not have {old_ig} in its neighborhood.")
        PP_dict[pp_id]["neighborhood"].discard(old_ig)
        pp_ig_edges[pp_id].discard(old_ig)
    for new_ig in new_igs:
        PP_dict[pp_id]["neighborhood"].add(new_ig)
        pp_ig_edges[pp_id].add(new_ig)

def return_boundary_IG(se_id, ig_id, SE_dict):
    # If the SE does not include ig_id in its nodes, that’s an error.
    nodes = list(SE_dict[se_id]["nodes"])
    if ig_id not in nodes:
        raise ValueError(f"return_boundary_IG: IG {ig_id} not found in SE {se_id} nodes {nodes}.")
    boundary_IG = nodes[0] if nodes[1] == ig_id else nodes[1]
    return boundary_IG



def perform_se_duplication(se_id, IG_dict, SE_dict, PSM_dict):
    """
    Duplicate an SE (edge) by inserting a new IG (internal node) at the midpoint
    between the two parent IGs. Only the local neighborhoods are updated:
      - The parent's neighborhood is updated by replacing the old SE ID with the new SE ID
        that connects the parent to the new IG.
      - The new IG’s neighborhood is initialized to contain just the two new SEs.
      - Any PSM that referenced the old SE in its FROM set is updated immediately.
    This approach avoids sweeping through all PP entries and keeps updates local.
    """
    global next_se_id  # using the global counter for new SE IDs

    if se_id not in SE_dict:
        raise KeyError(f"[ERROR] perform_se_duplication: SE {se_id} does not exist!")

    # Get the two parent IGs.
    parent_ig1, parent_ig2 = SE_dict[se_id]["nodes"]

    # Create the new IG at the midpoint of the two parents.
    new_ig_id = get_new_ig_id()
    new_pos = ((np.array(IG_dict[parent_ig1]["position"]) + 
                np.array(IG_dict[parent_ig2]["position"])) / 2).tolist()
    # Start with an empty neighborhood; we will add the new SE IDs next.
    IG_dict[new_ig_id] = {
        "position": new_pos,
        "neighborhood": set(),
        "psm_ids": set()
    }

    # Create two new SEs that connect each parent to the new IG.
    new_se1_id = get_new_se_id()
    new_se2_id = get_new_se_id()
    SE_dict[new_se1_id] = {
        "nodes": (parent_ig1, new_ig_id)
    }
    SE_dict[new_se2_id] = {
        "nodes": (new_ig_id, parent_ig2)
    }

    # Local update: For each parent IG, update its neighborhood.
    # Replace the reference to the old SE with the new SE that connects that parent to new_ig_id.
    update_neighbors(parent_ig1, {se_id}, {new_se1_id}, IG_dict)
    update_neighbors(parent_ig2, {se_id}, {new_se2_id}, IG_dict)

    # Set the new IG’s neighborhood to be exactly the two new SE IDs.
    IG_dict[new_ig_id]["neighborhood"].update({new_se1_id, new_se2_id})

    # Now update any PSM that was using se_id in its FROM set.
    # We assume that a PSM is located at one of the parent IGs (or along the SE),
    # so we replace the reference locally.
    for psm_key, psm in PSM_dict.items():
        if se_id in psm["FROM"]:
            current_psm_ig = psm["current_location"]
            if current_psm_ig == parent_ig1:
                replacement = new_se1_id
            elif current_psm_ig == parent_ig2:
                replacement = new_se2_id
            else:
                # If the PSM is not exactly on a parent, add both new SE IDs.
                replacement = None
            psm["FROM"].remove(se_id)
            if replacement is not None:
                psm["FROM"].add(replacement)
            else:
                psm["FROM"].update({new_se1_id, new_se2_id})
    # Finally, remove the original SE.
    del SE_dict[se_id]



def perform_se_reduction(se_id, IG_dict, SE_dict, PP_dict, PSM_dict, psm_queue):
    """
    Perform an SE reduction on the given SE id.
    (This is your existing SE reduction code.)
    At some point, if a reduction would cause the two IG endpoints to be the same,
    we call consolidate_multi_edges to remove the duplicate (self‐loop) edge.
    """
    # Example: get the two IG endpoints for the SE
    if se_id not in SE_dict:
        raise KeyError(f"[ERROR] SE {se_id} not found!")
    ig1, ig2 = SE_dict[se_id]["nodes"]
    if ig1 not in IG_dict or ig2 not in IG_dict:
        raise KeyError(f"[ERROR] One or both IGs {ig1}, {ig2} are missing!")

    # NEW: Bypass reduction if all PSMs on both IGs share se_id in their FROM set.
    if bypass_reduction_due_to_shared_edge(se_id, ig1, ig2, IG_dict, PSM_dict, SE_dict):
        return
    else:
        # For example, you would merge the two IGs, update PP connections, update PSMs, etc.
        #
        # At an appropriate point (after all the local updates and before returning),
        # you could call consolidate_multi_edges(IG_dict, SE_dict, PSM_dict)
        # to “clean up” any duplicate edges that may have been created during the reduction.
        #
        # For example:
        #
        #   ... perform local update operations ...
        #   consolidate_multi_edges(IG_dict, SE_dict, PSM_dict)
        #   print("[DEBUG] Completed SE reduction and consolidated multi-edges.")
        #
        # Then continue with any final steps or return.

        # (Rest of your reduction logic here.)

        # --- (1) Compute new IG’s position and initial neighborhood ---
        pos1 = np.array(IG_dict[ig1]["position"])
        pos2 = np.array(IG_dict[ig2]["position"])
        new_pos = ((pos1 + pos2) / 2).tolist()
        new_neighbors = (IG_dict[ig1]["neighborhood"] | IG_dict[ig2]["neighborhood"]) - {se_id}
        new_ig_id = get_new_ig_id()
        IG_dict[new_ig_id] = {"position": new_pos, "neighborhood": new_neighbors.copy(), "psm_ids": set()}

        # --- (2) Locally update PP and SE references from ig1 and ig2 ---
        for old_ig in (ig1, ig2):
            for elem in IG_dict[old_ig]["neighborhood"] - {se_id}:
                if elem in PP_dict:
                    update_pp_connection(elem, {old_ig}, {new_ig_id}, PP_dict, pp_ig_edges)
                if elem in SE_dict:
                    update_neighbors(elem, {old_ig}, {new_ig_id}, SE_dict)

        # --- (3) Partition the PSMs attached to ig1 and ig2 ---
        def partition_psms(ig):
            incoming = []
            outgoing = []
            for psm_id in IG_dict[ig]["psm_ids"]:
                if se_id in PSM_dict[psm_id]["FROM"]:
                    outgoing.append(psm_id)
                else:
                    incoming.append(psm_id)
            return incoming, outgoing

        incoming1, outgoing1 = partition_psms(ig1)
        incoming2, outgoing2 = partition_psms(ig2)

        def is_older(psm_a, psm_b):
            return int(psm_a[3:]) < int(psm_b[3:])

        psms_to_delete = set()

        # --- (4) Pair incoming from ig1 with outgoing from ig2 ---
        for psm_in in incoming1:
            if outgoing2:
                psm_out = outgoing2.pop(0)
                new_from = (PSM_dict[psm_in]["FROM"] | PSM_dict[psm_out]["FROM"]) - {se_id}
                PSM_dict[psm_out]["FROM"] = new_from
                psms_to_delete.add(psm_in)

        # --- (4b) Pair incoming from ig2 with outgoing from ig1 ---
        for psm_in in incoming2:
            if outgoing1:
                psm_out = outgoing1.pop(0)
                new_from = (PSM_dict[psm_in]["FROM"] | PSM_dict[psm_out]["FROM"]) - {se_id}
                PSM_dict[psm_out]["FROM"] = new_from
                psms_to_delete.add(psm_in)

        # --- (5) For any unpaired outgoing PSMs, update their FROM sets ---
        def get_bounding_element(ig, outgoing_from_set):
            candidates = IG_dict[ig]["neighborhood"] - {se_id}
            filtered = [cand for cand in candidates if cand not in outgoing_from_set]
            if len(filtered) == 1:
                return filtered[0]
            elif len(filtered) > 1:
                return filtered[0]
            else:
                return None

        for ig, psms in [(ig1, outgoing1), (ig2, outgoing2)]:
            for psm_id in psms:
                if psm_id in psms_to_delete:
                    continue
                partner_ig = ig2 if ig == ig1 else ig1
                bound = get_bounding_element(partner_ig, PSM_dict[psm_id]["FROM"])
                new_from = PSM_dict[psm_id]["FROM"].copy()
                new_from.discard(se_id)
                if bound is not None:
                    new_from.add(bound)
                PSM_dict[psm_id]["FROM"] = new_from
        # --- (6) Reassign all remaining (unmerged) PSMs to the new IG ---
        remaining_psm_ids = (IG_dict[ig1]["psm_ids"] | IG_dict[ig2]["psm_ids"]) - psms_to_delete
        for psm_id in remaining_psm_ids:
            PSM_dict[psm_id]["current_location"] = new_ig_id
            IG_dict[new_ig_id]["psm_ids"].add(psm_id)

        # --- (7) Remove merged (newer) PSMs ---
        for psm_id in psms_to_delete:
            if psm_id in psm_queue:
                psm_queue.remove(psm_id)
            if psm_id in PSM_dict:
                del PSM_dict[psm_id]

        # --- (8) Remove the old IGs and the reducing SE ---
        del IG_dict[ig1]
        del IG_dict[ig2]
        del SE_dict[se_id]

    

    return psm_queue


def perform_pp_split(PP_dict, IG_dict, SE_dict, PSM_dict, pp_to_split):
    """
    Split a passive partition (PP) into two new PPs by:
      1. Computing two new PP positions as offsets from the parent's position.
      2. Creating two new IGs at the same positions as the new PPs.
      3. Creating a new SE between these two new IGs.
      4. Updating each parent IG (i.e. every IG in the original PP’s neighborhood)
         by removing the old PP reference and adding the two new PP IDs.
      5. Recording new PP–IG edges locally.
      6. Creating new PSMs on each parent IG that remains connected.
      7. Finally, removing the original PP.
    
    This version only performs local (atomic) updates on the affected neighborhoods.
    """
    if pp_to_split not in PP_dict:
        return [], []
    
    # Get the parent's position and its directly connected IGs.
    parent_position = np.array(PP_dict[pp_to_split]["position"])
    parent_igs = list(PP_dict[pp_to_split]["neighborhood"])
    if not parent_igs:
        return [], []
    
    
    # Compute positions for the two new PP nodes.
    offset = np.array([1, 1, 0])
    child1_position = (parent_position - offset).tolist()
    child2_position = (parent_position + offset).tolist()
    
    # Create two new PP IDs.
    new_pp1_id = get_new_pp_id()
    new_pp2_id = get_new_pp_id()
    
    # Create two new IGs for the new PPs.
    new_ig1_id = get_new_ig_id()
    new_ig2_id = get_new_ig_id()
    
    IG_dict[new_ig1_id] = {"position": child1_position, "neighborhood": set(), "psm_ids": set()}
    IG_dict[new_ig2_id] = {"position": child2_position, "neighborhood": set(), "psm_ids": set()}
    
    # Create the new PP entries with neighborhoods that combine the parent's IGs and the new IG.
    PP_dict[new_pp1_id] = {
        "position": child1_position,
        "neighborhood": set(parent_igs) | {new_ig1_id},
    }
    PP_dict[new_pp2_id] = {
        "position": child2_position,
        "neighborhood": set(parent_igs) | {new_ig2_id},
    }
    
    # Update the passive edge mapping locally.
    pp_ig_edges[new_pp1_id] = set(parent_igs) | {new_ig1_id}
    pp_ig_edges[new_pp2_id] = set(parent_igs) | {new_ig2_id}
    
    # Update each parent's IG neighborhood: remove the old PP and add the two new PP IDs.
    for ig in parent_igs:
        if pp_to_split in IG_dict[ig]["neighborhood"]:
            IG_dict[ig]["neighborhood"].discard(pp_to_split)
        IG_dict[ig]["neighborhood"].update({new_pp1_id, new_pp2_id})
    
    # For each new IG, add its corresponding PP and also the new SE that will connect them.
    new_se_id = get_new_se_id()
    IG_dict[new_ig1_id]["neighborhood"].add(new_pp1_id)
    IG_dict[new_ig1_id]["neighborhood"].add(new_se_id)
    IG_dict[new_ig2_id]["neighborhood"].add(new_pp2_id)
    IG_dict[new_ig2_id]["neighborhood"].add(new_se_id)
    
    # Create the new SE connecting the two new IGs.
    SE_dict[new_se_id] = {
        "nodes": (new_ig1_id, new_ig2_id),
        "color": updated_edge_color
    }


    # For every parent IG, update or create a PSM to represent the new PP connection.
    new_psm_ids = []
    for ig in parent_igs:
        if ig in IG_dict:
            # The new FROM set for the split
            new_from_set = {new_pp1_id, new_pp2_id}
            # Get all propagating PSMs currently on the IG
            ig_psms = [psm_id for psm_id in IG_dict[ig]["psm_ids"]]
            # Find any that already reference the splitting PP (pp_to_split)
            psms_with_pp = [psm_id for psm_id in ig_psms if pp_to_split in PSM_dict[psm_id]["FROM"]]
            
            # If there are no propagating PSMs on this IG OR (if there is one and it does not contain pp_to_split)
            if not ig_psms or (len(ig_psms) == 1 and not psms_with_pp):
                # Create a new PSM for this IG.
                new_psm_id = get_new_psm_id()
                PSM_dict[new_psm_id] = {
                    "FROM": new_from_set.copy(),
                    "current_location": ig,
                    "propagating": True
                }
                new_psm_ids.append(new_psm_id)
                psm_queue.append(new_psm_id)
                IG_dict[ig]["psm_ids"].add(new_psm_id)
            else:
                # If at least one propagating PSM exists that already has the splitting PP in its FROM set,
                # update the first one we find.
                if psms_with_pp:
                    psm_to_update = psms_with_pp[0]
                    PSM_dict[psm_to_update]["FROM"].discard(pp_to_split)
                    PSM_dict[psm_to_update]["FROM"].update(new_from_set)
                else:
                    # Otherwise, create a new PSM.
                    new_psm_id = get_new_psm_id()
                    PSM_dict[new_psm_id] = {
                        "FROM": new_from_set.copy(),
                        "current_location": ig,
                        "propagating": True
                    }
                    new_psm_ids.append(new_psm_id)
                    psm_queue.append(new_psm_id)
                    IG_dict[ig]["psm_ids"].add(new_psm_id)

    # Finally, remove the old PP.
    del PP_dict[pp_to_split]
    
    # Return the new PSM IDs and the new SE ID for further propagation if needed.
    return new_psm_ids, [new_se_id]



# --- Helper Functions ---

def get_from_type(fs):
    if all(elem.startswith("p") for elem in fs):
        return "PP"
    elif all(elem.startswith("s") for elem in fs):
        return "SE"
    else:
        return "mixed"


def find_opposing_psm(current_ig, psm_id, IG_dict):
    for other in IG_dict[current_ig]["psm_ids"]:
        if other != psm_id:
            return other
    return None


def is_fused_ig(ig_id, IG_dict, PP_dict):
    """
    Returns True if all the neighbors of IG ig_id (as recorded in its 'neighborhood')
    are PP IDs.
    """
    neighbors = IG_dict[ig_id].get("neighborhood", set())
    if not neighbors:
        return False
    for elem in neighbors:
        if not elem.startswith("p"):
            return False
    return True

# --- Core Propagation Functions ---
def resolve_collision_psm(psm_prop, psm_opp, prop_psm_id, opp_psm_id, 
                          PSM_dict, IG_dict, SE_dict, PP_dict, psm_queue, next_psm_queue,
                          propagating_PSM_color, nodeOffset):
    """
    Resolve a collision between two PSMs using strictly local operations.
    
    The routine follows these steps:
      0. Preserve the current (central) IG as old_central.
      1. For each element in the propagating PSM’s FROM set, compute a new IG
         position with a radial offset and plan a local update (without applying it immediately).
         Immediately update local PP/SE connections from central IG to the new IG.
      2. For each element in the opposing PSM’s FROM set, if it is an SE, determine its
         bounding IG relative to old_central (using return_boundary_IG). If it is a PP,
         update its connection from the central IG to the new propagating IGs.
      3. For every bounding IG found, create new SE edges connecting it to each new propagating IG.
      4. Spawn new “opposing” PSMs on the new propagating IGs and new “propagating” PSMs on
         the bounding IGs.
      5. Remove the old opposing PSM and all local SE keys in its FROM set.
      6. Remove the central IG and its propagating PSM.
      7. Update local PSM queues with the newly spawned PSMs.
    
    All updates are performed locally so that the SERD model’s assumptions hold.
    """
    
    # (0) Preserve central IG.
    central_ig = psm_prop["current_location"]
    if central_ig not in IG_dict:
        return
    old_central = central_ig  # preserve original central IG for opposing SEs
    central_pos = np.array(IG_dict[central_ig]["position"])
    prop_from = sorted(list(psm_prop["FROM"]))
    opp_from  = sorted(list(psm_opp["FROM"]))
    opp_prop_bool = psm_opp["propagating"]
    # (1) For each element in the propagating FROM set, plan creation of a new IG.
    new_ig_by_elem = {}
    for idx, elem in enumerate(prop_from):
        angle = 2 * math.pi * idx / len(prop_from)
        offset = np.array([math.cos(angle), math.sin(angle), 0]) * nodeOffset
        new_ig_id = get_new_ig_id()
        new_ig_by_elem[elem] = new_ig_id
        # Immediately create the new IG locally.
        IG_dict[new_ig_id] = {
            "position": (central_pos + offset).tolist(),
            "neighborhood": {elem},
            "psm_ids": set()
        }
        # Update local connections from central_ig for this element.
        if elem in PP_dict:
            update_pp_connection(elem, {central_ig}, {new_ig_id}, PP_dict, pp_ig_edges)
        elif elem in SE_dict:
            update_neighbors(elem, {central_ig}, {new_ig_id}, SE_dict)
    new_ig_ids = list(new_ig_by_elem.values())

    # (2) Process opposing FROM set.
    bounding_IGs = set()
    opposing_PPs = set()
    for elem in opp_from:
        if elem.startswith("s") and elem in SE_dict:
            try:
                b_ig = return_boundary_IG(elem, old_central, SE_dict)
            except ValueError as e:
                b_ig = list(SE_dict[elem]["nodes"])[0]
            bounding_IGs.add(b_ig)
        elif elem.startswith("p") and elem in PP_dict:
            update_pp_connection(elem, {central_ig}, set(new_ig_ids), PP_dict, pp_ig_edges)
            opposing_PPs.add(elem)
            for ig in new_ig_ids:
                update_neighbors(ig, set(), {elem}, IG_dict)

    # (3) For each bounding IG, create new SE edges to each new propagating IG.
    new_se_edges = {}
    if bounding_IGs:
        for b_ig in bounding_IGs:
            for new_ig in new_ig_ids:
                new_se_id = get_new_se_id()
                new_se_edges[(b_ig, new_ig)] = new_se_id
                SE_dict[new_se_id] = {
                    "nodes": (new_ig, b_ig)
                }
                update_neighbors(new_ig, set(), {new_se_id}, IG_dict)
                update_neighbors(b_ig, set(), {new_se_id}, IG_dict)

    # (4) Spawn new opposing PSMs on each new propagating IG.
    new_opposing_psm_ids = {}
    for new_ig in new_ig_ids:
        local_SEs = {se for (b, n), se in new_se_edges.items() if n == new_ig} if bounding_IGs else set()
        new_from = local_SEs | opposing_PPs
        new_psm_id = get_new_psm_id()
        PSM_dict[new_psm_id] = {
            "FROM": new_from,
            "current_location": new_ig,
            "propagating": opp_prop_bool
        }
        IG_dict[new_ig]["psm_ids"].add(new_psm_id)
        new_opposing_psm_ids[new_ig] = new_psm_id
    
    # (5) Spawn new propagating PSMs on each bounding IG.
    new_propagating_psm_ids = {}
    if bounding_IGs:
        for b_ig in bounding_IGs:
            local_SEs = {se for (b, n), se in new_se_edges.items() if b == b_ig}
            new_psm_id = get_new_psm_id()
            PSM_dict[new_psm_id] = {
                "FROM": local_SEs,
                "current_location": b_ig,
                "propagating": True
            }
            IG_dict[b_ig]["psm_ids"].add(new_psm_id)
            new_propagating_psm_ids[b_ig] = new_psm_id

    # (6) Remove the old opposing PSM and any SE entries in its FROM set.
    if opp_psm_id in PSM_dict:
        del PSM_dict[opp_psm_id]
    for elem in opp_from:
        if elem.startswith("s"):
            del SE_dict[elem]
    
    # (7) Remove the central IG and its propagating PSM.
    if prop_psm_id in PSM_dict:
        del PSM_dict[prop_psm_id]
    if central_ig in IG_dict:
        del IG_dict[central_ig]
    
    # (8) Update local PSM queues.
    new_prop_set = set(new_propagating_psm_ids.values())
    if new_prop_set:
        next_psm_queue.extend(list(new_prop_set))
    new_opp_set = set(new_opposing_psm_ids.values())
    if opp_psm_id in psm_queue:
        idx = psm_queue.index(opp_psm_id)
        psm_queue[idx:idx+1] = list(new_opp_set)
    elif opp_psm_id in next_psm_queue:
        idx = next_psm_queue.index(opp_psm_id)
        next_psm_queue[idx:idx+1] = list(new_opp_set)
    else:
        next_psm_queue.extend(list(new_opp_set))
    return


def general_propagate_normal(psm, PSM_dict, IG_dict, SE_dict, PP_dict, next_psm_queue,
                             current_psm_id, nodeOffset):
    """
    Propagate a PSM from its current IG to a next local element using the local neighborhood.
    For each element in the FROM set, we now compute a unique offset based on that element’s type:
      - If the element is a PP, we use the PP’s position.
      - If the element is an SE, we get the bounding IG (the IG on the “other side” of the SE) and
        use its position.
    New IGs are created by offsetting the current IG’s position along these computed direction vectors.
    """
    current_ig = psm["current_location"]
    from_set = psm["FROM"]
    local_neigh = IG_dict[current_ig].get("neighborhood", set())
    
    # Step 1: Choose a candidate neighbor not already in the FROM set.
    candidate = None
    candidate_type = None
    for nb in local_neigh:
        if nb not in from_set:
            if nb in SE_dict:
                candidate = nb
                candidate_type = "SE"
                break
            elif nb in PP_dict:
                candidate = nb
                candidate_type = "PP"
                break
    if candidate is None:
        IG_dict[current_ig]["psm_ids"].discard(current_psm_id)
        if current_psm_id in PSM_dict:
            del PSM_dict[current_psm_id]
        return

    # Step 2: If candidate is a PP, perform local PP partitioning.
    if candidate_type == "PP":
        current_pos = np.array(IG_dict[current_ig]["position"])
        new_ig_ids = set()
        for idx, elem in enumerate(from_set):
            angle = 2 * math.pi * idx / len(from_set)
            offset = np.array([math.cos(angle), math.sin(angle), 0]) * nodeOffset
            new_ig_id = get_new_ig_id()
            new_ig_ids.add(new_ig_id)
            IG_dict[new_ig_id] = {
                "position": (current_pos + offset).tolist(),
                "neighborhood": {elem, candidate},
                "psm_ids": set()
            }
            if elem.startswith("p") and elem in PP_dict:
                if current_ig in PP_dict[elem]["neighborhood"]:
                    update_pp_connection(elem, {current_ig}, {new_ig_id}, PP_dict, pp_ig_edges)
            elif elem.startswith("s") and elem in SE_dict:
                if current_ig in SE_dict[elem]["nodes"]:
                    update_neighbors(elem, {current_ig}, {new_ig_id}, SE_dict)
        update_pp_connection(candidate, {current_ig}, new_ig_ids, PP_dict, pp_ig_edges)
        del IG_dict[current_ig]
        if current_psm_id in PSM_dict:
            del PSM_dict[current_psm_id]
        return

    # Step 3: Candidate is a SE.
    se = SE_dict[candidate]
    nodes = se["nodes"]
    if nodes[0] == current_ig:
        next_elem = nodes[1]
    elif nodes[1] == current_ig:
        next_elem = nodes[0]
    else:
        return

    # Step 4: For each element in the FROM set, determine its target position.
    # Instead of using the same target for all, we choose:
    #   - For PP elements: use PP's position.
    #   - For SE elements: use the position of the bounding IG (the IG opposite to current_ig).
    current_pos = np.array(IG_dict[current_ig]["position"])
    new_ig_ids = {}
    for elem in from_set:
        if elem in PP_dict:
            target_pos = np.array(PP_dict[elem]["position"])
        elif elem in SE_dict:
            # Use the bounding IG of this SE relative to current_ig.
            try:
                bound_ig = return_boundary_IG(elem, current_ig, SE_dict)
                target_pos = np.array(IG_dict[bound_ig]["position"])
            except Exception as e:
                target_pos = np.array(IG_dict[next_elem]["position"])
        else:
            target_pos = np.array(IG_dict[next_elem]["position"])

        diff = target_pos - current_pos
        norm = np.linalg.norm(diff)
        if norm == 0:
            direction_vector = np.array([nodeOffset, 0, 0])
        else:
            direction_vector = diff / norm
        offset = direction_vector * nodeOffset
        new_ig_id = get_new_ig_id()
        new_ig_ids[elem] = new_ig_id
        IG_dict[new_ig_id] = {
            "position": (current_pos + offset).tolist(),
            "neighborhood": set(),  # to be updated below
            "psm_ids": set()
        }

    # Step 5: Update local connections for each new IG.
    for elem, new_ig in new_ig_ids.items():
        if elem in PP_dict:
            update_pp_connection(elem, {current_ig}, {new_ig}, PP_dict, pp_ig_edges)
            IG_dict[new_ig]["neighborhood"].add(elem)
        elif elem in SE_dict:
            update_neighbors(elem, {current_ig}, {new_ig}, SE_dict)
            IG_dict[new_ig]["neighborhood"].add(elem)
    
    # Step 6: Create new SE edges from each new IG to next_elem.
    new_se_ids = set()
    for new_ig in new_ig_ids.values():
        new_se_id = get_new_se_id()
        new_se_ids.add(new_se_id)
        SE_dict[new_se_id] = {
            "nodes": (new_ig, next_elem)
        }
        IG_dict[new_ig]["neighborhood"].add(new_se_id)
    update_neighbors(next_elem, {candidate}, new_se_ids, IG_dict)
    if candidate in SE_dict:
        del SE_dict[candidate]
    
    # Remove the old IG.
    del IG_dict[current_ig]
    if next_elem in IG_dict:
        IG_dict[next_elem]["psm_ids"].add(current_psm_id)
    psm["current_location"] = next_elem
    psm["FROM"] = new_se_ids
    next_psm_queue.append(current_psm_id)
    return


def propagate_psm_main(psm_id, PSM_dict, IG_dict, SE_dict, PP_dict, 
                         psm_queue, next_psm_queue, frame_counter, propagating_PSM_color, nodeOffset):
    if psm_id not in PSM_dict:
        return
    psm = PSM_dict[psm_id]
    current_ig = psm["current_location"]
    if current_ig not in IG_dict:
        return

    if len(IG_dict[current_ig]["psm_ids"]) > 1:
        opposing_psm_id = find_opposing_psm(current_ig, psm_id, IG_dict)
        if opposing_psm_id is None:
            return
        opposing_psm = PSM_dict[opposing_psm_id]
        resolve_collision_psm(psm, opposing_psm, psm_id, opposing_psm_id,
                              PSM_dict, IG_dict, SE_dict, PP_dict,
                              psm_queue, next_psm_queue, propagating_PSM_color, nodeOffset)
    else:
        general_propagate_normal(psm, PSM_dict, IG_dict, SE_dict, PP_dict,
                                 next_psm_queue, psm_id, nodeOffset)


def propagate_psm(PSM_dict, IG_dict, SE_dict, PP_dict, current_psm_id, psm_queue, next_psm_queue, frame_counter, propagating_PSM_color, nodeOffset):
    propagate_psm_main(current_psm_id, PSM_dict, IG_dict, SE_dict, PP_dict, psm_queue, next_psm_queue, frame_counter, propagating_PSM_color, nodeOffset)





# Define fission operation ----------
def perform_pp_fission_on_IG(ig_id, PP_dict, IG_dict, SE_dict, PSM_dict, pp_ig_edges):
    """
    Perform PP fission on an IG that is connected to exactly two PPs.
    That is, if IG ig_id has exactly two PP neighbors, then remove ig_id
    and replace it with two new IGs (each connected to one PP) joined by a new SE.
    
    Updates:
      - Creates two new IGs (n_new1 and n_new2) with positions slightly offset from the original.
      - Creates a new SE connecting these two IGs.
      - In each PP (from the original PP neighbors) replaces the reference to ig_id with the appropriate new IG.
      - Reassigns any PSMs that were attached to ig_id to one of the new IGs (using a simple heuristic).
      - Deletes the old IG.
    """
    # Retrieve the IG data.
    ig = IG_dict.get(ig_id)
    if ig is None:
        print(f"[ERROR] IG {ig_id} not found for PP fission.")
        return

    # Identify PP neighbors: those elements in ig["neighborhood"] that are PP IDs.
    pp_neighbors = {elem for elem in ig["neighborhood"] if elem in PP_dict}
    if len(pp_neighbors) != 2:
        print(f"[WARNING] IG {ig_id} does not have exactly 2 PP neighbors (found {len(pp_neighbors)}); fission skipped.")
        return

    pp_list = list(pp_neighbors)
    p1, p2 = pp_list[0], pp_list[1]

    # Get the original IG position.
    orig_pos = np.array(ig["position"])
    # For a sensible offset, we can use the vector between the two PPs.
    p1_pos = np.array(PP_dict[p1]["position"])
    p2_pos = np.array(PP_dict[p2]["position"])
    # Compute a direction from p1 to p2; if they coincide, choose an arbitrary direction.
    dir_vec = p2_pos - p1_pos
    if np.linalg.norm(dir_vec) == 0:
        dir_vec = np.array([1, 0, 0])
    else:
        dir_vec = dir_vec / np.linalg.norm(dir_vec)
    # Compute a perpendicular vector for splitting.
    arbitrary = np.array([0, 0, 1])
    perp = np.cross(dir_vec, arbitrary)
    if np.linalg.norm(perp) == 0:
        perp = np.array([0, 1, 0])
    else:
        perp = perp / np.linalg.norm(perp)
    # Use a small offset for the new positions.
    offset = 0.1 * perp  
    new_pos1 = (orig_pos + offset).tolist()
    new_pos2 = (orig_pos - offset).tolist()

    # Create two new IGs.
    new_ig1 = get_new_ig_id()
    new_ig2 = get_new_ig_id()
    IG_dict[new_ig1] = {"position": new_pos1, "neighborhood": {p1}, "psm_ids": set()}
    IG_dict[new_ig2] = {"position": new_pos2, "neighborhood": {p2}, "psm_ids": set()}

    # Create a new SE connecting the two new IGs.
    new_se = get_new_se_id()
    SE_dict[new_se] = {"nodes": (new_ig1, new_ig2), "color": updated_edge_color, "last_updated": 0}
    IG_dict[new_ig1]["neighborhood"].add(new_se)
    IG_dict[new_ig2]["neighborhood"].add(new_se)

    # Update the PP connections.
    # For PP p1, replace ig_id with new_ig1.
    if ig_id in PP_dict[p1]["neighborhood"]:
        PP_dict[p1]["neighborhood"].remove(ig_id)
        PP_dict[p1]["neighborhood"].add(new_ig1)
    if p1 in pp_ig_edges and ig_id in pp_ig_edges[p1]:
        pp_ig_edges[p1].remove(ig_id)
        pp_ig_edges[p1].add(new_ig1)
    # For PP p2, replace ig_id with new_ig2.
    if ig_id in PP_dict[p2]["neighborhood"]:
        PP_dict[p2]["neighborhood"].remove(ig_id)
        PP_dict[p2]["neighborhood"].add(new_ig2)
    if p2 in pp_ig_edges and ig_id in pp_ig_edges[p2]:
        pp_ig_edges[p2].remove(ig_id)
        pp_ig_edges[p2].add(new_ig2)

    # Reassign any PSMs that were attached to the old IG.
    for psm_id in ig.get("psm_ids", set()):
        psm = PSM_dict.get(psm_id)
        if psm is None:
            continue
        # As a simple heuristic, if the PSM's FROM set contains p2, assign to new_ig2; else assign to new_ig1.
        if p2 in psm.get("FROM", set()):
            target = new_ig2
        else:
            target = new_ig1
        psm["current_location"] = target
        IG_dict[target]["psm_ids"].add(psm_id)

    # Finally, remove the original IG.
    del IG_dict[ig_id]


def perform_pp_fission_update(update_el_id, PP_dict, IG_dict, SE_dict, PSM_dict, pp_ig_edges):
    """
    Check if the update element is an IG eligible for PP fission (i.e. has exactly two PP neighbors).
    If so, perform PP fission on it.
    """
    if update_el_id in IG_dict:
        ig = IG_dict[update_el_id]
        # Determine the set of PP neighbors.
        pp_neighbors = {elem for elem in ig.get("neighborhood", set()) if elem in PP_dict}
        if len(pp_neighbors) == 2:
            perform_pp_fission_on_IG(update_el_id, PP_dict, IG_dict, SE_dict, PSM_dict, pp_ig_edges)




# ─── 1) single-trajectory “latent” runner ──────────────────────────────────
def simulate_latent(T_latent,
                    SE0, IG0, PP0,
                    pp_split_p0, dup_p, red_p,
                    decay_gamma, t_stop_splits,
                    embed_freq=0.0):
    """
    Single‐trajectory latent run:
      - Phase 1 (per TS): non‐deterministic SE dup/reduction + PP splits
      - Phase 2 (per TS): deterministic, one‐hop PSM propagation
    """
    import copy, random, math

    # --- 1) deep-copy your starting network
    SE  = copy.deepcopy(SE0)
    IG  = copy.deepcopy(IG0)
    PP  = copy.deepcopy(PP0)

    # --- 2) seed exactly one PSM per PP→IG boundary
    PSM = {}
    psm_queue = []
    def _new_psm_id():
        return f"psm{len(PSM)}"
    for pp_id, pp in PP.items():
        for ig_id in pp["neighborhood"]:
            pid = _new_psm_id()
            PSM[pid] = {
                "FROM": {pp_id},
                "current_location": ig_id,
                "propagating": True
            }
            IG[ig_id]["psm_ids"].add(pid)
            psm_queue.append(pid)

    next_psm_queue = []

    # small helper for optionally re-embedding
    def _maybe_embed():
        if embed_freq>0 and random.random()<embed_freq:
            spring_electric_embedding_local(IG, SE, PP)

    # --- 3) run T_latent steps, each with two phases
    for t in range(T_latent):
        # ─── PHASE 1 ─── schedule & apply exactly one SE‐dup/reduction per SE
        # build & shuffle update queue: all idle‐SEs + all PPs
        uq = [eid for eid,se in SE.items() if 'pending_op' not in se] + list(PP.keys())
        random.shuffle(uq)

        # current split-prob
        p_s = pp_split_p0 * math.exp(-decay_gamma * t)
        if t > t_stop_splits:
            p_s = 0.0

        # 3.1) decide & mark ops
        for eid in uq:
            if eid in SE:
                op = random.choices(['dup','red','none'], [dup_p, red_p, 1-(dup_p+red_p)])[0]
                if op!='none':
                    SE[eid]['pending_op'] = op
                _maybe_embed()
            else:  # eid is a PP
                if random.random() < p_s:
                    perform_pp_split(PP, IG, SE, PSM, eid)
                _maybe_embed()

        # 3.2) apply _all_ pending SE ops immediately
        for sid in list(SE):
            if 'pending_op' in SE[sid]:
                if SE[sid]['pending_op']=='dup':
                    perform_se_duplication(sid, IG, SE, PSM)
                else:
                    perform_se_reduction(sid, IG, SE, PP, PSM, psm_queue)
                del SE[sid]['pending_op']
                _maybe_embed()

        # ─── PHASE 2 ─── one synchronous sweep of _all_ PSMs
        next_psm_queue.clear()
        for pid in psm_queue:
            propagate_psm(PSM, IG, SE, PP,
                          pid, psm_queue, next_psm_queue,
                          t, propagating_PSM_color, nodeOffset)
            _maybe_embed()
        # swap into next round
        psm_queue, next_psm_queue = next_psm_queue, []

    return SE, IG, PP, PSM


# ─── 2) MC-ensemble + 3D+2D histories ─────────────────────────────────────
def run_ensemble(M, T_max,
                 SE0, IG0, PP0, PSM0,
                 pp_split_p0, dup_p, red_p,
                 decay_gamma, t_stop_splits,
                 embed_freq, grid_res):
    # history containers
    pp_history   = [set() for _ in range(T_max)]
    ig_history   = [set() for _ in range(T_max)]
    edge_history = [set() for _ in range(T_max)]

    # scalar metrics (12)
    metrics = {
        'pp_count':   np.zeros((T_max, M)),
        'ig_count':   np.zeros((T_max, M)),
        'se_count':   np.zeros((T_max, M)),
        'avg_radius': np.zeros((T_max, M)),
        'density':    np.zeros((T_max, M)),
        'p_split':    np.zeros((T_max, M)),
        'sep':        np.zeros((T_max, M)),
        'ig_deg':     np.zeros((T_max, M)),
        'pp_deg':     np.zeros((T_max, M)),
        'topo_dim':   np.zeros((T_max, M)),
        'clust':      np.zeros((T_max, M)),
        'psm_count':  np.zeros((T_max, M)),
    }

    # 2D-grid for density fields
    xedges    = np.linspace(-axisScale, axisScale, grid_res+1)
    yedges    = np.linspace(-axisScale, axisScale, grid_res+1)
    pp_hist2d = np.zeros((T_max, grid_res, grid_res), dtype=float)
    ig_hist2d = np.zeros((T_max, grid_res, grid_res), dtype=float)

    for run in tqdm(range(M), desc="MC runs"):
        SE  = copy.deepcopy(SE0)
        IG  = copy.deepcopy(IG0)
        PP  = copy.deepcopy(PP0)
        PSM = copy.deepcopy(PSM0)

        # seed your PSM queue once per run
        psm_queue      = list(PSM.keys())
        next_psm_queue = []

        for t in trange(T_max, desc="  steps", leave=False):
            # 0) compute current split-prob
            p_s = pp_split_p0 * math.exp(-decay_gamma * t)
            if t > t_stop_splits:
                p_s = 0.0
            metrics['p_split'][t, run] = p_s

            # ─── PHASE 1: non-deterministic updates ────────────────────────────────────
            uq = [eid for eid,se in SE.items() if 'pending_op' not in se] + list(PP.keys())
            random.shuffle(uq)

            for eid in uq:
                if eid in SE and 'pending_op' not in SE[eid]:
                    op = random.choices(['dup','red','none'], [dup_p, red_p, 1-(dup_p+red_p)])[0]
                    if op!='none':
                        SE[eid]['pending_op'] = op
                        SE[eid]['last_updated'] = t
                    if random.random() < embed_freq:
                        spring_electric_embedding_local(IG, SE, PP)

                elif eid in PP and random.random() < p_s:
                    new_psms, _ = perform_pp_split(PP, IG, SE, PSM, eid)
                    # inject them immediately into the next propagation phase
                    psm_queue.extend(new_psms)
                    if random.random() < embed_freq:
                        spring_electric_embedding_local(IG, SE, PP)

            # apply *all* pending SE ops right away
            for sid, se in list(SE.items()):
                if 'pending_op' in se:
                    if se['pending_op']=='dup':
                        perform_se_duplication(sid, IG, SE, PSM)
                    else:
                        perform_se_reduction(sid, IG, SE, PP, PSM, psm_queue)
                    if random.random() < embed_freq:
                        spring_electric_embedding_local(IG, SE, PP)

            # ─── PHASE 2: deterministic PSM propagation ─────────────────────────────────
            next_psm_queue.clear()

            # guard against infinite loops
            max_iters = max(100, 10 * len(PSM))
            iters = 0

            # keep going until nothing left to propagate this TS
            while psm_queue:
                iters += 1
                if iters > max_iters:
                    print(f"[WARN] PSM propagation exceeded {max_iters} iterations at TS={t}")
                    break

                pid = psm_queue.pop(0)
                propagate_psm(PSM, IG, SE, PP,
                              pid, psm_queue, next_psm_queue,
                              t, propagating_PSM_color, nodeOffset)

                if random.random() < embed_freq:
                    spring_electric_embedding_local(IG, SE, PP)

            # swap into the next TS
            psm_queue, next_psm_queue = next_psm_queue, []

            # ─── rest of your per-TS bookkeeping (embed, histograms, metrics…) ─────────
            if random.random() < embed_freq:
                spring_electric_embedding_local(IG, SE, PP)

            # 6) positions & edges
            pp_pts = np.array([PP[k]['position'] for k in PP])
            ig_pts = np.array([IG[k]['position'] for k in IG])

            # cloud separation via KMeans
            if pp_pts.shape[0] >= 2:
                km = KMeans(n_clusters=2, random_state=0).fit(pp_pts)
                c1, c2 = km.cluster_centers_
                sep = np.linalg.norm(c1 - c2)
            else:
                sep = 0.0
            metrics['sep'][t, run] = sep

            for p in pp_pts: pp_history[t].add(tuple(p))
            for p in ig_pts: ig_history[t].add(tuple(p))

            # record edges
            for sid, se in SE.items():
                n1, n2 = se['nodes']
                p1, p2 = IG[n1]['position'], IG[n2]['position']
                edge_history[t].add((tuple(p1), tuple(p2)))
            for pp_id, pp in PP.items():
                ppos = tuple(pp['position'])
                for ig_id in pp['neighborhood']:
                    edge_history[t].add((ppos, tuple(IG[ig_id]['position'])))

            # 6b) 2D histograms
            if pp_pts.size:
                Hpp,_,_ = np.histogram2d(pp_pts[:,0], pp_pts[:,1],
                                         bins=[xedges,yedges])
                pp_hist2d[t] += Hpp
            if ig_pts.size:
                Hig,_,_ = np.histogram2d(ig_pts[:,0], ig_pts[:,1],
                                         bins=[xedges,yedges])
                ig_hist2d[t] += Hig

            # 7) scalar metrics
            ppN, igN, seN, psmN = len(PP), len(IG), len(SE),len(psm_queue)
            metrics['pp_count'][t,run]  = ppN
            metrics['ig_count'][t,run]  = igN
            metrics['se_count'][t,run]  = seN
            metrics['psm_count'][t,run] = psmN

            centroid = pp_pts.mean(axis=0) if pp_pts.size else np.zeros(3)
            dists    = np.linalg.norm(pp_pts - centroid, axis=1) \
                          if pp_pts.size else np.array([0.])
            R_rms    = np.sqrt((dists**2).mean()) if dists.size else 0.
            R_max    = dists.max() if dists.size else 1.0
            metrics['avg_radius'][t,run] = R_rms
            vol = (4/3)*math.pi*(R_max**3) if R_max>0 else 1.0
            metrics['density'][t,run]    = ppN/vol

            # avg degrees
            metrics['ig_deg'][t,run] = np.mean([len(IG[n]['neighborhood']) for n in IG]) 
                                       
            metrics['pp_deg'][t,run] = np.mean([len(PP[p]['neighborhood']) for p in PP]) 
                                       

            # fractal dimension
            D,_,_ = estimate_topological_dimension(IG, SE)
            metrics['topo_dim'][t,run] = D

            # build the full SERD graph (IG + PP nodes, SE + PP-IG edges)
            G = nx.Graph()
            G.add_nodes_from(IG.keys())
            G.add_nodes_from(PP.keys())
            # add the SE edges between IGs
            G.add_edges_from(se["nodes"] for se in SE.values())
            # add the PP→IG edges
            for pp_id, pp in PP.items():
                for ig_id in pp["neighborhood"]:
                    G.add_edge(pp_id, ig_id)

    return pp_history, ig_history, edge_history, metrics, (xedges, yedges, pp_hist2d, ig_hist2d)


# ─── 3) pick latent vs straight-chain and run ensemble ──────────────────


if use_latent:
    SE_init, IG_init, PP_init, PSM_init = simulate_latent(
        latent_steps,
        SE_dict, IG_dict, PP_dict, PSM_dict,
        pp_split_probability,
        duplication_probability,
        reduction_probability,
        decay_gamma,
        t_stop_splits,
        embed_freq
    )
else:
    # Kick‐off your chain‐of‐10 at start of script
    initialize_chain(n=init_SE_filament)
    SE_init, IG_init, PP_init, PSM_init = SE_dict, IG_dict, PP_dict, PSM_dict

# run MC ensemble
pp_history, ig_history, edge_history, mets, (xedges, yedges, pp2d, ig2d) = run_ensemble(
    M, T_max,
    SE_init, IG_init, PP_init, PSM_init,
    pp_split_probability,
    duplication_probability,
    reduction_probability,
    decay_gamma,
    t_stop_splits,
    embed_freq,
    grid_res=50
)


# ─── Precompute means & IQRs ──────────────────────────────────────────────
ts = np.arange(T_max)

def median_q25_q75(arr):
    m  = np.median(arr, axis=1)
    lo = np.percentile(arr, 25, axis=1)
    hi = np.percentile(arr, 75, axis=1)
    return m, lo, hi

pp_mean = mets['pp_count'].mean(axis=1)
pp_std  = mets['pp_count'].std(axis=1)

ig_m, ig_lo, ig_hi    = median_q25_q75(mets['ig_count'])
se_m, se_lo, se_hi    = median_q25_q75(mets['se_count'])
r_m,  r_lo,  r_hi     = median_q25_q75(mets['avg_radius'])
d_m,  d_lo,  d_hi     = median_q25_q75(mets['density'])
ps_m, ps_lo, ps_hi    = median_q25_q75(mets['p_split'])
sep_m, sep_lo, sep_hi = median_q25_q75(mets['sep'])
igdeg_m, igdeg_lo, igdeg_hi = median_q25_q75(mets['ig_deg'])
ppdeg_m, ppdeg_lo, ppdeg_hi = median_q25_q75(mets['pp_deg'])
topo_m, topo_lo, topo_hi    = median_q25_q75(mets['topo_dim'])
clust_m, clust_lo, clust_hi = median_q25_q75(mets['clust'])
psm_mean = mets['psm_count'].mean(axis=1)


# ─── 4) BUILD 2×2 FIGURE ─────────────────────────────────────────────────
fig = plt.figure(figsize=(18,12))
outer = GridSpec(2, 2, figure=fig, wspace=0.3, hspace=0.3)

# 4a) 3D sim (top-left)
ax_sim = fig.add_subplot(outer[0,0], projection='3d')
scatter_pp = ax_sim.scatter([],[],[],
    c='black', s=pp_size, alpha=pp_alpha, label='PP')
scatter_ig = ax_sim.scatter([],[],[],
    c='orange', s=ig_size, alpha=ig_alpha, label='IG')
edge_coll = Line3DCollection([], colors='black',
    linewidths=edge_width, alpha=edge_alpha)
ax_sim.add_collection(edge_coll)
ax_sim.set_xlim(-axisScale,axisScale)
ax_sim.set_ylim(-axisScale,axisScale)
ax_sim.set_zlim(-axisScale,axisScale)
ax_sim.set_title(f"SERD Simulation (M={M})")
ax_sim.legend(loc='upper left')

# 4b) small histograms (bottom-left)
hist_gs = outer[1,0].subgridspec(2,1, hspace=0.4)
ax_hpp = fig.add_subplot(hist_gs[0,0])
pcm_pp = ax_hpp.pcolormesh(xedges, yedges, pp2d[0].T,
                           cmap='Blues', vmin=0, vmax=pp2d.max(),
                           shading='auto')
ax_hpp.set_title("⟨PP density⟩"); ax_hpp.set_xlabel("x"); ax_hpp.set_ylabel("y")
ax_hpp.scatter([],[], c='red', s=80, marker='x', label='centroids')
ax_hpp.legend()

ax_hig = fig.add_subplot(hist_gs[1,0])
pcm_ig = ax_hig.pcolormesh(xedges, yedges, ig2d[0].T,
                           cmap='Oranges', vmin=0, vmax=ig2d.max(),
                           shading='auto')
ax_hig.set_title("⟨IG density⟩"); ax_hig.set_xlabel("x"); ax_hig.set_ylabel("y")
ax_hig.scatter([],[], c='red', s=80, marker='x', label='centroids')
ax_hig.legend()

# lock both histogram axes to your global axisScale
for ax in (ax_hpp, ax_hig):
    ax.set_xlim(-axisScale, axisScale)
    ax.set_ylim(-axisScale/2, axisScale/2)
    ax.set_aspect('equal', 'box')   # optional, makes x/y units equal

# 4c) time-series metrics (right half split into 2×(3×2))
metrics_outer = outer[:,1].subgridspec(2,1, hspace=0.5)
# top-right 6 panels
top_gs = metrics_outer[0].subgridspec(3,2, hspace=0.8, wspace=0.5)
axs_top = [fig.add_subplot(top_gs[i,j]) for i in range(3) for j in range(2)]
# bottom-right 6 panels
bot_gs = metrics_outer[1].subgridspec(3,2, hspace=0.8, wspace=0.5)
axs_bot = [fig.add_subplot(bot_gs[i,j]) for i in range(3) for j in range(2)]
axs_ts = axs_top + axs_bot

labels = [
    "# PP", "# IG", "# SE",
    "RMS radius", "density", "p_split",
    "cloud sep.", "avg IG deg", "avg PP deg",
    "D_topo", "avg clust", "PSM count"
]
series = [
    pp_mean,    ig_m,    se_m,
    r_m,        d_m,     ps_m,
    sep_m,      igdeg_m, ppdeg_m,
    topo_m,     clust_m, psm_mean
]
iqrs = [
    (pp_mean-pp_std, pp_mean+pp_std),
    (ig_lo, ig_hi),   (se_lo, se_hi),
    (r_lo, r_hi),     (d_lo, d_hi),
    (ps_lo, ps_hi),   (sep_lo, sep_hi),
    (igdeg_lo, igdeg_hi), (ppdeg_lo, ppdeg_hi),
    (topo_lo, topo_hi),   (clust_lo, clust_hi),
    None
]

lines = []
for ax, lbl, data, iqr in zip(axs_ts, labels, series, iqrs):
    ax.set_title(lbl)
    ln, = ax.plot([], [], 'k-')
    if iqr is not None:
        lo, hi = iqr
        ax.fill_between(ts, lo, hi, color='gray', alpha=0.3)
    ax.set_xlim(0, T_max-1)
    lines.append(ln)

axs_top[0].set_ylabel("value")


# ─── 5) update function & animation ──────────────────────────────────────
def update_frame(t):
    # 3D points & edges
    pts_pp = np.array(list(pp_history[t]))
    scatter_pp._offsets3d = ((pts_pp[:,0], pts_pp[:,1], pts_pp[:,2]) 
                              if pts_pp.size else ([],[],[]))
    pts_ig = np.array(list(ig_history[t]))
    scatter_ig._offsets3d = ((pts_ig[:,0], pts_ig[:,1], pts_ig[:,2]) 
                              if pts_ig.size else ([],[],[]))
    segs = [[p1,p2] for (p1,p2) in edge_history[t]]
    edge_coll.set_segments(segs)

    # time-series
    for ln, data in zip(lines, series):
        ln.set_data(ts[:t+1], data[:t+1])

    # histograms
    pcm_pp.set_array(pp2d[t].T.ravel())
    pcm_ig.set_array(ig2d[t].T.ravel())

    # centroids on histograms
    if pts_pp.shape[0] >= 2:
        cen_pp = KMeans(n_clusters=2, random_state=0).fit(pts_pp[:,:2]).cluster_centers_
    elif pts_pp.shape[0]==1:
        cen_pp = pts_pp[:,:2]
    else:
        cen_pp = np.empty((0,2))
    ax_hpp.collections[-1].set_offsets(cen_pp)

    if pts_ig.shape[0] >= 2:
        cen_ig = KMeans(n_clusters=2, random_state=0).fit(pts_ig[:,:2]).cluster_centers_
    elif pts_ig.shape[0]==1:
        cen_ig = pts_ig[:,:2]
    else:
        cen_ig = np.empty((0,2))
    ax_hig.collections[-1].set_offsets(cen_ig)

    return [scatter_pp, scatter_ig, edge_coll] + lines + [pcm_pp, pcm_ig]

ani = FuncAnimation(fig, update_frame, frames=T_max, interval=200, blit=False)
plt.show()