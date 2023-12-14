
#include "global.h"
#include "helper.h"
#include "index_base.h"
#include "masstree/masstree.hh"
#include "masstree/kvthread.hh"
#include "masstree/masstree_tcursor.hh"
#include "masstree/masstree_insert.hh"
#include "masstree/masstree_get.hh"
#include "masstree/masstree_split.hh"
#include "masstree/masstree_remove.hh"
#include "masstree/masstree_scan.hh"
#include "masstree/str.hh"




typedef struct mass_node {
    // TODO bad hack!
    void **pointers; // for non-leaf nodes, point to mass_nodes
    bool is_leaf;
    idx_key_t *keys;
    mass_node *parent;
    UInt32 num_keys;
    mass_node *next;
    bool latch;
    pthread_mutex_t locked;
    latch_t latch_type;
    UInt32 share_cnt;
    std::atomic<uint64_t> t_access;
} mass_node;

struct glob_param {
	uint64_t part_id;
};

class index_masstree : public index_base {
public:
	RC			init(uint64_t part_cnt);
	RC			init(uint64_t part_cnt, table_t * table);
	bool 		index_exist(idx_key_t key); // check if the key exist. 
	RC 			index_insert(idx_key_t key, itemid_t * item, int part_id = -1);
	RC	 		index_read(idx_key_t key, itemid_t * &item, 
					uint64_t thd_id, int64_t part_id = -1);
	RC	 		index_read(idx_key_t key, itemid_t * &item, int part_id = -1);
	RC	 		index_read(idx_key_t key, itemid_t * &item);
	RC 			index_next(uint64_t thd_id, itemid_t * &item, bool samekey = false);

private:
	// index structures may have part_cnt = 1 or PART_CNT.
	uint64_t part_cnt;
	RC			make_lf(uint64_t part_id, mass_node *& node);
	RC			make_nl(uint64_t part_id, mass_node *& node);
	RC		 	make_node(uint64_t part_id, mass_node *& node);
	
	RC 			start_new_tree(glob_param params, idx_key_t key, itemid_t * item);
	RC 			find_leaf(glob_param params, idx_key_t key, idx_acc_t access_type, mass_node *& leaf, mass_node  *& last_ex);
	RC 			find_leaf(glob_param params, idx_key_t key, idx_acc_t access_type, mass_node *& leaf);
	RC			insert_into_leaf(glob_param params, mass_node * leaf, idx_key_t key, itemid_t * item);
	// handle split
	RC 			split_lf_insert(glob_param params, mass_node * leaf, idx_key_t key, itemid_t * item);
	RC 			split_nl_insert(glob_param params, mass_node * node, UInt32 left_index, idx_key_t key, mass_node * right);
	RC 			insert_into_parent(glob_param params, mass_node * left, idx_key_t key, mass_node * right);
	RC 			insert_into_new_root(glob_param params, mass_node * left, idx_key_t key, mass_node * right);

	int			leaf_has_key(mass_node * leaf, idx_key_t key);
	
	UInt32 		cut(UInt32 length);
	UInt32	 	order; // # of keys in a node(for both leaf and non-leaf)
	mass_node ** 	roots; // each partition has a different root
	mass_node *   find_root(uint64_t part_id);

	bool 		latch_node(mass_node * node, latch_t latch_type);
	latch_t		release_latch(mass_node * node);
	RC		 	upgrade_latch(mass_node * node);
	// clean up all the LATCH_EX up tp last_ex
	RC 			cleanup(mass_node * node, mass_node * last_ex);

	// the leaf and the idx within the leaf that the thread last accessed.
	mass_node *** cur_leaf_per_thd;
	UInt32 ** 		cur_idx_per_thd;
};

