function result = hlp_memavail()
% Get the amount of potentially available memory, in bytes
% This is usually more than what is reported by hlp_memfree, because some memory is tentatively 
% allocated by the OS for caches, etc.
bean = java.lang.management.ManagementFactory.getOperatingSystemMXBean();
result = bean.getTotalPhysicalMemorySize - bean.getCommittedVirtualMemorySize;
