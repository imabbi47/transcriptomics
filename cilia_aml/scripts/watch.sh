#!/usr/bin/env bash
# Emits ONE event when the pipeline reaches a terminal state, its process
# disappears (teardown), so silence never masks a dead job.
BASE=/home/abhi/github-projects/brainstroming_transcriptomics/cilia_aml/data
LOG=$BASE/pipeline.log
while true; do
  if grep -qE "PIPELINE DONE|FATAL" "$LOG" 2>/dev/null; then
    echo "TERMINAL: $(grep -E 'PIPELINE DONE|FATAL' "$LOG" | tail -1)"
    break
  fi
  # [p] trick avoids the grep matching its own command line
  if ! ps -eo args 2>/dev/null | grep -q "[p]ipeline.sh"; then
    echo "PIPELINE PROCESS GONE without a DONE marker — likely torn down; re-launch needed"
    break
  fi
  sleep 120
done
echo "==== last 30 lines of pipeline.log ===="
tail -30 "$LOG"
