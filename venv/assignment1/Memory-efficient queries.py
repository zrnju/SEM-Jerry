import cellxgene_census
import tiledbsoma

with cellxgene_census.open_soma() as census:

    human = census["census_data"]["homo_sapiens"]
    query = human.axis_query(
       measurement_name = "RNA",
       obs_query = tiledbsoma.AxisQuery(
           value_filter = "tissue == 'brain' and sex == 'male'"
       )
    )

    # Continued below
    # Continued from above

    iterator = query.X("raw").tables()

    # Get an iterative slice as pyarrow.Table
    raw_slice = next (iterator)
    ...

    # Continued from above
    query.close()