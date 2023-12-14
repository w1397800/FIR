FIR
-----------------

FIR is a  a NVM OLTP engine with fast rebuilding of the DRAM index, achieving instant system recovery while maintaining high runtime throughput.




Experiment Environment
------------

Experiments are conducted on a server equipped with an Intel Xeon Gold 5218R CPU (20 physical cores @2.10GHz and 27.5 MB LLC) and 96 GB DDR4 DRAM, running Ubuntu 20.04.3 LTS. Each core supports two hardware threads, resulting in a total of 40 threads. There are 768 GB (128x6 GB) Intel Optane DC Persistent Memory NVDIMMs in the system. Optane PM has two accessibility modes, Memory Mode and APP Direct mode. We chose APP Direct mode for easy mapping of DRAM and NVM into the same virtual address space. We deploy file systems in fs-dax mode on NVM, followed by employing libpmem within PMDK to establish mappings between NVM files and a process's virtual memory. To ensure data persistence to NVM, we utilize \textit{clwb} and \textit{sfence}. The entire codebase is developed in C/C++ and compiled using gcc 7.5.0.


Experimental Results
-------------

- **Transaction Performance**
    - TPCC Performance: https://github.com/w1397800/FIR/wiki/TPCC-Performance-with-Different-OLTP-Engines
    - YCSB Performance:  
- **Recovery Performance:**
    - Recovery with Unplanned Crash: 
    - Recovery with Manual Shutdown: 
- **Index Checkpoint Performance:**
    - TPCC Performance with Different OLTP Engines: 
    - Throughput over 100 Seconds: 
