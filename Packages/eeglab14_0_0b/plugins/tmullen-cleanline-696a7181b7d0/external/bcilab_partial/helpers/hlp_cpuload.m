function result = hlp_cpuload
% Get the relative CPU load
osbean = java.lang.management.ManagementFactory.getOperatingSystemMXBean();
result = osbean.getSystemLoadAverage() / osbean.getAvailableProcessors();
