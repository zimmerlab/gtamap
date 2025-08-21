package server

import "github.com/KleinSamuel/gtamap/src/config"

func GetGenomicIntervalFromRecord(r *EnhancedRecord) *Interval {

	lastRegion, err := r.MappedGenome.GetLastRegion()
	if !err {
		return nil
	}

	return &Interval{
		Contig: r.Rname,
		Start:  r.Pos,
		End:    lastRegion.End,
	}
}

// IsInTargetRegion checks if the given record is within the target region
// defined in the configuration.
// 0 = no overlap
// 1 = contained entirely
// 2 = overlaps only left
// 3 = overlaps only right
// 4 = target is contained in record (overlaps right and left)
func IsInTargetRegion(r *EnhancedRecord) int {

	contig := r.Rname

	if contig != config.GetTargetContig() {
		return 0
	}

	contigStart := r.Pos

	lastRegion, err := r.MappedGenome.GetLastRegion()
	if !err {
		return 0
	}
	contigEnd := lastRegion.End

	if contigStart < config.GetTargetStart() && contigEnd > config.GetTargetEnd() {
		return 4 // target is contained in record
	}
	if contigStart >= config.GetTargetStart() && contigEnd <= config.GetTargetEnd() {
		return 1 // record is contained entirely in target
	}
	if contigStart < config.GetTargetStart() && contigEnd > config.GetTargetStart() {
		return 2 // record overlaps only left
	}
	if contigStart < config.GetTargetEnd() && contigEnd > config.GetTargetEnd() {
		return 3 // record overlaps only right
	}

	return -1
}

func ContainedInTargetContig(contig string, start int, end int) bool {
	if contig != config.GetTargetContig() {
		return false
	}
	if start < config.GetTargetStart() || end > config.GetTargetEnd() {
		return false
	}
	return true
}
