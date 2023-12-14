/* Masstree
 * Eddie Kohler, Yandong Mao, Robert Morris
 * Copyright (c) 2012-2014 President and Fellows of Harvard College
 * Copyright (c) 2012-2014 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, subject to the conditions
 * listed in the Masstree LICENSE file. These conditions include: you must
 * preserve this copyright notice, and you cannot mention the copyright
 * holders in advertising related to the Software without their permission.
 * The Software is provided WITHOUT ANY WARRANTY, EXPRESS OR IMPLIED. This
 * notice is a summary of the Masstree LICENSE file; the license in that file
 * is legally binding.
 */
#ifndef MASSTREE_INSERT_HH
#define MASSTREE_INSERT_HH
#include "masstree_get.hh"
#include "masstree_split.hh"
#include "storage/peloton/store/storage_manager.h"

namespace Masstree {

template <typename P>
bool tcursor<P>::find_insert(threadinfo& ti)
{
    find_locked(ti);
    original_n_ = n_;
    original_v_ = n_->full_unlocked_version_value();

//    if(kx_.p==9){
//
//      int a = 0;
//      a++;
//    }
//    if(kx_.p>=14){
//
//      int a = 0;
//      a++;
//    }
//    if(kx_.p==17){
//
//      int a = 0;
//      a++;
//    }
//    if(ka_.ikey() == 7719732711266451456){
//      int a = 0;
//      a++;
//    }

    // maybe we found it
    if (state_)
        return true;

    // otherwise mark as inserted but not present
    state_ = 2;

    // maybe we need a new layer
    if (kx_.p >= 0){
//      StorageManager::GetInstance()->addCountLayer();
      return make_new_layer(ti);
    }

    // mark insertion if we are changing modification state
    if (unlikely(n_->modstate_ != leaf<P>::modstate_insert)) {
        masstree_invariant(n_->modstate_ == leaf<P>::modstate_remove);
        n_->mark_insert();
        n_->modstate_ = leaf<P>::modstate_insert;
    }

    // try inserting into this node
    if (n_->size() < n_->width) {
        kx_.p = permuter_type(n_->permutation_).back();
        // don't inappropriately reuse position 0, which holds the ikey_bound
        if (likely(kx_.p != 0) || !n_->prev_ || n_->ikey_bound() == ka_.ikey()) {
            n_->assign(kx_.p, ka_, ti);
//            if(kx_.p != kx_.i){
//              int a = 0;
//              a++;
//            }
//            if(kx_.p ==0&& kx_.i ==0 ){
//              int a = 0;
//              a++;
//            }
            return false;
        }
    }

    // otherwise must split
    return make_split(ti);
}

//
//template <typename P>
//node_base<P> *tcursor<P>::add_leaf(node_base<P>* parent, threadinfo& ti){
//  n_ = (static_cast<const leaf<P> *>(parent));
//  original_n_ = n_;
//  original_v_ = n_->full_unlocked_version_value();
//
//  // inserting into this node
//  if (n_->size() < n_->width) {
//    kx_.p = permuter_type(n_->permutation_).back();
//      n_->assign(kx_.p, ka_, ti);
//      return false;
//  }
//
//  return make_split(ti);
//}

template <typename P>
void tcursor<P>::add_son(node_base<P>* parent, threadinfo& ti, int p, int mid, int height, key_type ka){

//  leaf_type* n_;
//  key_type ka_;
//  key_indexed_position kx_;
//  node_base<P>* root_;

  auto *iparent =static_cast<leaf_type*>(parent);

  n_ = iparent;
  ka_ = ka;
  root_ = parent;
  kx_.p = p;
  kx_.i = p;

  n_->assign(p, ka, ti);
  finish_insert();

}

template <typename P>
void tcursor<P>::add_son(leaf_type* iparent, threadinfo& ti, int p,
                         int mid, int height, key_type ka){

  //  leaf_type* n_;
  //  key_type ka_;
  //  key_indexed_position kx_;
  //  node_base<P>* root_;

//  auto *iparent =static_cast<leaf_type*>(parent);

  n_ = iparent;
  ka_ = ka;
  root_ = iparent;
  kx_.p = p;
  kx_.i = p;

  n_->assign(p, ka, ti);
  finish_insert();

}
//
//template <typename P>
//void tcursor<P>::add_layer(node_base<P>* parent, threadinfo& ti, int p,
//    int mid, int height, std::vector<chpt_tuple> chpt_arr,  std::vector<uint8_t> pre_array, key_type ka){
//
//  auto *iparent =static_cast<leaf_type*>(parent);
//
//  leaf_type* nl = leaf_type::make_root(0, iparent, ti);
//
//  iparent->lv_[kx_.p] = nl;
//
//
//
//}

//template <typename P>
//node_base<P> *tcursor<P>::build_tree(node_base<P>* parent, threadinfo& ti, uint start, uint end, uint height){
//  if (start > end) return nullptr;
//
//  for (int i = 1; i <= 15; ++i) {
//    int mid = (start + end) / 2 + (i - 8);
//
//
//    if(height == 0){
//      add_leaf(parent, ti);
//    }else {
//      int p = i-1;
//      add_node(parent, ti);
//    }
//  }
//
//
//  TreeNode* root = new TreeNode((start + end) / 2, height);
//  for (int i = 1; i <= 15; ++i) {
//    int mid = (start + end) / 2 + (i - 8);
//
//
//
//    root->children.push_back(build_tree(mid + 1, end, height + 1));
//  }
//
//}


template <typename P>
bool tcursor<P>::make_new_layer(threadinfo& ti) {
    key_type oka(n_->ksuf(kx_.p));
    ka_.shift();
    int kcmp = oka.compare(ka_);

    // Create a twig of nodes until the suffixes diverge
    leaf_type* twig_head = n_;
    leaf_type* twig_tail = n_;
    while (kcmp == 0) {
        leaf_type* nl = leaf_type::make_root(0, twig_tail, ti);
        nl->assign_initialize_for_layer(0, oka);
        if (twig_head != n_)
            twig_tail->lv_[0] = nl;
        else
            twig_head = nl;
        nl->permutation_ = permuter_type::make_sorted(1);
        twig_tail = nl;
        new_nodes_.emplace_back(nl, nl->full_unlocked_version_value());
        oka.shift();
        ka_.shift();
        kcmp = oka.compare(ka_);
    }

    // Estimate how much space will be required for keysuffixes
    size_t ksufsize;
    if (ka_.has_suffix() || oka.has_suffix())
        ksufsize = (std::max(0, ka_.suffix_length())
                    + std::max(0, oka.suffix_length())) * (n_->width / 2)
            + n_->iksuf_[0].overhead(n_->width);
    else
        ksufsize = 0;
    leaf_type *nl = leaf_type::make_root(ksufsize, twig_tail, ti);
    nl->assign_initialize(0, kcmp < 0 ? oka : ka_, ti);
    nl->assign_initialize(1, kcmp < 0 ? ka_ : oka, ti);
    nl->lv_[kcmp > 0] = n_->lv_[kx_.p];
    nl->lock(*nl, ti.lock_fence(tc_leaf_lock));
    if (kcmp < 0)
        nl->permutation_ = permuter_type::make_sorted(1);
    else {
        permuter_type permnl = permuter_type::make_sorted(2);
        permnl.remove_to_back(0);
        nl->permutation_ = permnl.value();
    }
    // In a prior version, recursive tree levels and true values were
    // differentiated by a bit in the leafvalue. But this constrains the
    // values users could assign for true values. So now we use bits in
    // the key length, and changing a leafvalue from true value to
    // recursive tree requires two writes. How to make this work in the
    // face of concurrent lockless readers? Mark insertion so they
    // retry.
    n_->mark_insert();
    fence();
    if (twig_tail != n_)
        twig_tail->lv_[0] = nl;
    if (twig_head != n_)
        n_->lv_[kx_.p] = twig_head;
    else
        n_->lv_[kx_.p] = nl;
    n_->keylenx_[kx_.p] = n_->layer_keylenx;
    updated_v_ = n_->full_unlocked_version_value();
    n_->unlock();
    n_ = nl;
    kx_.i = kx_.p = kcmp < 0;
    return false;
}

template <typename P>
void tcursor<P>::finish_insert()
{
    permuter_type perm(n_->permutation_);
    masstree_invariant(perm.back() == kx_.p);
    perm.insert_from_back(kx_.i);
    fence();
    n_->permutation_ = perm.value();
}

template <typename P>
inline void tcursor<P>::finish(int state, threadinfo& ti)
{
    if (state < 0 && state_ == 1) {
        if (finish_remove(ti))
            return;
    } else if (state > 0 && state_ == 2)
        finish_insert();
    // we finally know this!
    if (n_ == original_n_)
        updated_v_ = n_->full_unlocked_version_value();
    else
        new_nodes_.emplace_back(n_, n_->full_unlocked_version_value());
    n_->unlock();
}

} // namespace Masstree
#endif
