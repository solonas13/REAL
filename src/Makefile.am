noinst_HEADERS = \
	acgtnMap.hpp \
	AsynchronousBufferReader.hpp \
	AsynchronousReader.hpp \
	AsynchronousWriter.hpp \
	AutoArray.hpp \
	AutoTextArray.hpp \
	BinRadixSort.hpp \
	CharBuffer.hpp \
	CharTermTable.hpp \
	ComputeScore.hpp \
	countReads.hpp \
	ERank222B.hpp \
	FastAReader.hpp \
	FastDecoder.hpp \
	FastDecodeInfo.hpp \
	FastFileDecoder.hpp \
	FastFileDecoderBase.hpp \
	FastFileDecoderState.hpp \
	FastIDBlock.hpp \
	FastIDDecoder.hpp \
	FastIDFileDecoder.hpp \
	FASTQEntry.hpp \
	FastQReader.hpp \
	FastQualityDecoder.hpp \
	FastQualityFileDecoder.hpp \
	FastQualityIDFileDecoder.hpp \
	FastQualitySubDecoder.hpp \
	FastReaderBase.hpp \
	FastSubDecoder.hpp \
	gappedMatching.hpp \
	getFileID.hpp \
	getFileList.hpp \
	getHistSize.hpp \
	getLookupTable.hpp \
	getNumLists.hpp \
	getPatternSig.hpp \
	getPhysicalMemory.hpp \
	getRadixSortTemp.hpp \
	getSampleBits.hpp \
	getText.hpp \
	IPow.hpp \
	ListSet.hpp \
	ListSetBlockReader.hpp \
	Lock.hpp \
	MappedMagic.hpp \
	MapTextFile.hpp \
	Mask.hpp \
	match.hpp \
	matchAll.hpp \
	matchAllImplementation.cpp \
	MatcherBase.hpp \
	matchUnique.hpp \
	matchUniqueImplementation.cpp \
	ParallelRadixSort.hpp \
	parseLines.hpp \
	Pattern.hpp \
	PatternBaseBlock.hpp \
	PatternBlock.hpp \
	PatternIdReader.hpp \
	PatternQualityBaseBlock.hpp \
	PopCountTable.hpp \
	popcnt.hpp \
	PosixMutex.hpp \
	PosixSemaphore.hpp \
	radixsort.hpp \
	RangeSet.hpp \
	RangeVector.hpp \
	RealOptions.hpp \
	RealTimeClock.hpp \
	ReorderFastA.hpp \
	ReorderFastQ.hpp \
	RestMatch.hpp \
	RestWordBuffer.hpp \
	Scoring.hpp \
	SignatureConstruction.hpp \
	SpaceTable.hpp \
	StreamTextFile.hpp \
	stringFunctions.hpp \
	SynchronousQueue.hpp \
	SynchronousReader.hpp \
	TemporaryFile.hpp \
	TerminatableSynchronousQueue.hpp \
	types.hpp \
	u4_sort.hpp \
	u8_sort.hpp \
	UniqueMatchInfo.hpp \
	u_sort.hpp

noinst_LIBRARIES = libscoring.a

noinst_PROGRAMS=genpat randstr reorderPat printScores
bin_PROGRAMS=real

libscoring_a_SOURCES = Scoring.cpp
libscoring_a_CXXFLAGS = @FFLOATSTORE@

singleinst=matchAllFastA32nonsse4noscores.cpp matchAllFastA32nonsse4scores.cpp matchAllFastA32sse4noscores.cpp matchAllFastA32sse4scores.cpp matchAllFastA64nonsse4noscores.cpp matchAllFastA64nonsse4scores.cpp matchAllFastA64sse4noscores.cpp matchAllFastA64sse4scores.cpp matchAllFastQ32nonsse4noscores.cpp matchAllFastQ32nonsse4scores.cpp matchAllFastQ32sse4noscores.cpp matchAllFastQ32sse4scores.cpp matchAllFastQ64nonsse4noscores.cpp matchAllFastQ64nonsse4scores.cpp matchAllFastQ64sse4noscores.cpp matchAllFastQ64sse4scores.cpp matchUniqueFastA32nonsse4noscores.cpp matchUniqueFastA32nonsse4scores.cpp matchUniqueFastA32sse4noscores.cpp matchUniqueFastA32sse4scores.cpp matchUniqueFastA64nonsse4noscores.cpp matchUniqueFastA64nonsse4scores.cpp matchUniqueFastA64sse4noscores.cpp matchUniqueFastA64sse4scores.cpp matchUniqueFastQ32nonsse4noscores.cpp matchUniqueFastQ32nonsse4scores.cpp matchUniqueFastQ32sse4noscores.cpp matchUniqueFastQ32sse4scores.cpp matchUniqueFastQ64nonsse4noscores.cpp matchUniqueFastQ64nonsse4scores.cpp matchUniqueFastQ64sse4noscores.cpp matchUniqueFastQ64sse4scores.cpp
	
real_SOURCES=\
	countReads.cpp \
	FastDecodeInfo.cpp \
	getFileID.cpp \
	getFileList.cpp \
	getNumLists.cpp \
	getPhysicalMemory.cpp \
	PopCountTable.cpp \
	real.cpp \
	RealOptions.cpp \
	StaticInitialization.cpp \
	${singleinst}

real_LDADD = libscoring.a

genpat_SOURCES=countReads.cpp genpat.cpp PopCountTable.cpp StaticInitialization.cpp
randstr_SOURCES=randstr.cpp
reorderPat_SOURCES=reorderPat.cpp RealOptions.cpp getPhysicalMemory.cpp FastDecodeInfo.cpp
reorderPat_LDADD = libscoring.a
printScores_SOURCES=printScores.cpp RealOptions.cpp getPhysicalMemory.cpp
printScores_LDADD = libscoring.a
