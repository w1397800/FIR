FIR
-----------------

FIR is an NVM OLTP engine with fast rebuilding of the DRAM index, achieving instant system recovery while maintaining high runtime throughput.




Experiment Environment
------------

Experiments were conducted on a server equipped with an Intel Xeon Gold 5218R CPU (20 physical cores @2.10GHz and 27.5 MB LLC) and 256 GB DDR4 DRAM, running Ubuntu 20.04.3 LTS. Each core supports two hardware threads, resulting in a total of 40 threads. There are 768 GB (128x6 GB) Intel Optane DC Persistent Memory NVDIMMs in the system. Optane PM has two accessibility modes, Memory Mode and APP Direct mode. We choose APP Direct mode for easy mapping of DRAM and NVM into the same virtual address space. We deploy file systems in fs-dax mode on NVM, followed by employing libpmem within PMDK to establish mappings between NVM files and a process's virtual memory. To ensure data persistence to NVM, we utilize clwb and sfence. The entire codebase is developed in C/C++ and compiled using gcc 7.5.0.

Experimental Results
-------------

- **Transaction Performance**
    - TPCC Performance: https://github.com/w1397800/FIR/wiki/TPCC-Performance-with-Different-OLTP-Engines
    - YCSB Performance: https://github.com/w1397800/FIR/wiki/YCSB-Performance
- **Recovery Performance:**
    - Recovery with Unplanned Crash: https://github.com/w1397800/FIR/wiki/Recovery-with-Unplanned-Crash
    - Recovery with Manual Shutdown: https://github.com/w1397800/FIR/wiki/Recovery-with-Manual-Shutdown
- **Index Checkpoint Performance:**
    - Checkpoint Speed & Checkpoint Scale: https://github.com/w1397800/FIR/wiki/Checkpoint-Speed-&-Checkpoint-Scale
    - Throughput over 100 Seconds: https://github.com/w1397800/FIR/wiki/Throughput-over-100-Seconds
