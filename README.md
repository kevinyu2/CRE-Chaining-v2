# CRE-Chaining-v2

## Chaining

### chaining.py

Contains driver functions for chaining (considering both possible orientations of the sequence): 
- chain_driver(anchors, is_weighted) -> global chaining, takes in list of [(a, c)] or [(a, c, weight)]
- chain_driver_np(anchors_np, is_weighted) -> takes in a np array with shape (n, 2) or (n, 3) for weighted
- chain_local_driver(anchors, match, mismatch, gap, is_weighted) -> if is_weighted, match score ignored
- chain_local_driver_np(anchors_np, match, mismatch, gap, is_weighted) -> if is_weighted, match score ignored

Contains individual chaining functions (only condiering given orientation)
- chain(anchors)
- chain_np(anchors_np)
- chain_local(anchors, match, mismatch, gap)
- chain_local_np(anchors_np, match, mismatch, gap)
- chain_weighted(anchors)
- chain_weighted_np(anchors_np)
- chain_local_weighted(anchors, match, mismatch, gap)
- chain_local_weighted_np(anchors_np, match, mismatch, gap)

Detailed info in the file header