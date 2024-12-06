import os
import json
import sqlite3


def load_prs_metadata(json_path):
    """
    Loads the prs_meta.json file.

    :param json_path: Path to the JSON file.
    :return: Dictionary containing PRS metadata.
    """
    try:
        # Read and parse the JSON file
        with open(json_path, "r") as file:
            return json.load(file)

    except Exception as e:
        print(f"Error loading prs_meta.json: {e}")
        return {}


def get_beta_stats(db_path, table_name):
    """
    Calculates the min and max beta values, the maximum possible score, and the minimum possible score.

    :param db_path: Path to the SQLite database.
    :param table_name: Name of the table to query.
    :return: Tuple (min_beta, max_beta, max_score, min_score).
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # Fetch min and max beta
        query_min_max = f"SELECT MIN(beta), MAX(beta) FROM {table_name}"
        cursor.execute(query_min_max)
        min_beta, max_beta = cursor.fetchone()

        # Calculate max possible score: sum(beta * 2) where beta > 0
        query_max_score = f"SELECT SUM(beta * 2) FROM {table_name} WHERE beta > 0"
        cursor.execute(query_max_score)
        max_score = cursor.fetchone()[0] or 0

        # Calculate min possible score: sum(beta * 2) where beta < 0
        query_min_score = f"SELECT SUM(beta * 2) FROM {table_name} WHERE beta < 0"
        cursor.execute(query_min_score)
        min_score = cursor.fetchone()[0] or 0

        return min_beta, max_beta, max_score, min_score

    except sqlite3.Error as e:
        print(f"Error accessing database table '{table_name}': {e}")
        return None, None, None, None

    finally:
        if conn:
            conn.close()


def get_hla_int_stats(db_path, table_name):
    """
    Calculates the overall minimum and maximum scores for the db_int table.

    :param db_path: Path to the SQLite database.
    :param table_name: Name of the db_int table to query.
    :return: Tuple (overall_min, overall_max).
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # Case 1: Rows where both a1 and a2 are NOT NULL
        query_not_null = f"""
            SELECT a1, a2, beta
            FROM {table_name}
            WHERE a1 IS NOT NULL AND a2 IS NOT NULL AND a1 != 'NULL' AND a2 != 'NULL'
        """
        cursor.execute(query_not_null)
        rows_not_null = cursor.fetchall()

        # Extract pairs from Non-NULL rows
        non_null_pairs = {(row[0], row[1]) for row in rows_not_null}

        # Case 2: Rows with NULL values representing duplicates or unique pairs
        query_null = f"""
            SELECT COALESCE(a1, a2) AS resolved_a, COALESCE(a2, a1) AS resolved_b, beta * 2
            FROM {table_name}
            WHERE a1 IS NULL OR a2 IS NULL OR a1 = 'NULL' OR a2 = 'NULL'
        """
        cursor.execute(query_null)
        rows_null = cursor.fetchall()

        # Filter rows_null to exclude duplicates already present in
        # non_null_pairs
        unique_null_rows = [
            row for row in rows_null
            # Resolve (NULL, DQ81) to (DQ81, DQ81)
            if (row[1], row[1]) not in non_null_pairs
        ]

        # Find the min and max rows from both sets
        all_rows = []
        for row in rows_not_null:
            # beta as the key, full row as value
            all_rows.append((row[2], row))
        for row in unique_null_rows:
            # beta * 2 as the key, full row as value
            all_rows.append((row[2], row))

        if all_rows:
            overall_min = min(row[0] for row in all_rows)
            overall_max = max(row[0] for row in all_rows)
        else:
            overall_min, overall_max = None, None

        return overall_min, overall_max

    except sqlite3.Error as e:
        print(f"Error accessing database table '{table_name}': {e}")
        return None, None

    finally:
        if conn:
            conn.close()


def update_json_with_min_max(json_path, updated_data):
    """
    Updates the JSON file with min and max values for each score.

    :param json_path: Path to the JSON file.
    :param updated_data: Dictionary with updated min and max values for each score.
    """
    try:
        # Load existing JSON data
        with open(json_path, "r") as file:
            data = json.load(file)

        # Update each score with its min and max values
        for score_name, values in updated_data.items():
            if score_name in data:
                data[score_name]["min"] = values.get("min")
                data[score_name]["max"] = values.get("max")
            else:
                print(f"Warning: {score_name} not found in the JSON file.")

        # Save the updated data back to the file
        with open(json_path, "w") as file:
            json.dump(data, file, indent=4)

        print(f"Updated JSON file saved at {json_path}.")

    except Exception as e:
        print(f"Error updating JSON file: {e}")


def analyze_and_update_json(json_path, db_path):
    """
    Analyze PRS scores, calculate min/max for each score, and update the JSON file.

    :param json_path: Path to the JSON file.
    :param db_path: Path to the SQLite database.
    """
    # Load PRS metadata
    prs_metadata = load_prs_metadata(json_path)
    if not prs_metadata:
        print("No metadata loaded. Exiting.")
        return

    updated_data = {}
    for score_name, meta in prs_metadata.items():
        method = meta.get("method", "unknown")
        db_table = meta.get("db_table")

        if method in ["additive", "grouped"] and db_table:
            # Fetch beta statistics
            min_beta, max_beta, max_score, min_score = get_beta_stats(
                db_path, db_table)
            updated_data[score_name] = {
                "min": min_score,
                "max": max_score,
            }

        elif method == "hla_int" and db_table and "db_int" in meta:
            # Process db_table as additive
            min_main, max_main, main_max_score, main_min_score = get_beta_stats(
                db_path, db_table)

            # Process db_int for HLA interactions
            min_hla, max_hla = get_hla_int_stats(db_path, meta["db_int"])

            # Calculate overall scores
            if main_min_score is not None and main_max_score is not None and min_hla is not None and max_hla is not None:
                overall_min = main_min_score + min_hla
                overall_max = main_max_score + max_hla
                updated_data[score_name] = {
                    "min": overall_min,
                    "max": overall_max,
                }

        else:
            print(f"Skipping {score_name} (Method: {method})")

    # Update the JSON file with the calculated min/max values
    update_json_with_min_max(json_path, updated_data)


if __name__ == "__main__":
    # Locate the JSON and database files
    script_location = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(script_location, "prs_meta.json")
    db_path = os.path.join(script_location, "../", "SQL", "variants.db")

    # Analyze PRS scores and update JSON
    analyze_and_update_json(json_path, db_path)
