From FASTQ to matrix
====================


Mapping
-------


.. currentmodule:: pytadbit.mapping.full_mapper

.. autofunction:: full_mapping


Quality check and plotting
--------------------------

.. currentmodule:: pytadbit.mapping.analyze

.. autofunction:: hic_map

.. autofunction:: insert_sizes

.. autofunction:: correlate_matrices

.. autofunction:: plot_strand_bias_by_distance

.. autofunction:: eig_correlate_matrices

.. autofunction:: plot_distance_vs_interactions

.. autofunction:: plot_iterative_mapping

.. autofunction:: plot_genomic_distribution


.. currentmodule:: pytadbit.utils.fastq_utils

.. autofunction:: quality_plot


Filtering
---------

.. currentmodule:: pytadbit.mapping.restriction_enzymes

.. autofunction:: map_re_sites

.. autofunction:: repaired
   
.. currentmodule:: pytadbit.mapping

.. autofunction:: get_intersection

.. autofunction:: merge_2d_beds

.. currentmodule:: pytadbit.mapping.filter

.. autofunction:: filter_reads

.. autofunction:: apply_filter

