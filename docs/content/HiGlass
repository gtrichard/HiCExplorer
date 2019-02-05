Visualizing HiCExplorer data in HiGlass
=======================================

Exporting h5 to cool
--------------------

It's possible to export Hi-C Matrices produced by HiCExplorer to
`bioconductor <http://bioconductor.org/>`__ in R, which allows us to use
existing bioconductor infrastructure for differential Hi-C analysis. The
tool **hicConvertFormat** allows us to write Hi-C matrices in a format that can
easily be imported in bioconductor as **GInteractions** object. Below
is an example.

```
hicConvertFormat
```

Balancing the cool matrix
-------------------------

```
cooler balance
```

Creating a multi resolution balanced cool matrix (mcool)
--------------------------------------------------------

```
cooler zoomify
```

Installing and launching HiGlass
--------------------------------


```
install higlass-manage
install docker
higlass-manage start
```

Visualizing the mcool matrix in HiGlass
---------------------------------------

```
higlass-manage ingest
higlass-manage view
```
