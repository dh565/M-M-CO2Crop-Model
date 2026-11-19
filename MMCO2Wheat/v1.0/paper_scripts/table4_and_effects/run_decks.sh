#!/bin/bash
# run_decks.sh <exe> <decks dir> <out dir> [workers]: run every deck of <decks dir> into <out dir>/<deck>/{input,output,run.log}
EXE=$1; DECKS=$2; OUT=$3; W=${4:-8}
mkdir -p "$OUT"
run_one() {
  d=$1; w="$OUT/$d"; rm -rf "$w"; mkdir -p "$w/output"; cp -R "$DECKS/$d/input" "$w/input"
  ( cd "$w" && "$EXE" > run.log 2>&1; echo "exit=$?" >> run.log )
  echo "$d $(tail -1 "$w/run.log")"
}
export -f run_one; export EXE DECKS OUT
ls "$DECKS" | grep -v -E "\.(csv|txt)$" | xargs -P "$W" -I{} bash -c 'run_one {}'
echo "RUNS DONE $(date)"
