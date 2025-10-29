def format_time(elapsed_seconds):
  # Calculate hours, minutes, and remaining seconds
  hours = int(elapsed_seconds // 3600)
  minutes = int((elapsed_seconds % 3600) // 60)
  seconds = int(elapsed_seconds % 60)

  # Format the output
  return f"{hours:02d}:{minutes:02d}:{seconds:02d}"
